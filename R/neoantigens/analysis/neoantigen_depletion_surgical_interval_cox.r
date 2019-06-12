#Survival analyses
#-------------------------

library(odbc)
library(DBI)
library(ggplot2)
library(survival)

con <- DBI::dbConnect(odbc::odbc(), "VerhaakDB")  

q <- "WITH roman_to_int (grade, grade_int)
AS
(
     VALUES ('I', 1),
            ('II', 2),
            ('III', 3),
            ('IV', 4)
)
SELECT
    tp.tumor_pair_barcode,
    tp.case_barcode,
    tp.tumor_barcode_a,
    tp.tumor_barcode_b,
	clin.idh_codel_subtype AS subtype,
	cas1.case_age_diagnosis_years AS age,
    (CASE
     WHEN s1.surgery_location = s2.surgery_location AND (s1.surgery_laterality = s2.surgery_laterality OR (s1.surgery_laterality IS NULL AND s2.surgery_laterality IS NULL)) THEN 'Local'
     WHEN s1.surgery_location <> s2.surgery_location OR s1.surgery_laterality <> s2.surgery_laterality THEN 'Distal'
     END) AS recurrence_location,
    (SELECT bool_or(treatment_tmz) FROM clinical.surgeries ss WHERE ss.case_barcode = tp.case_barcode AND ss.surgery_number >= s1.surgery_number AND ss.surgery_number < s2.surgery_number) AS received_tmz,
    (SELECT sum(treatment_tmz_cycles) FROM clinical.surgeries ss WHERE ss.case_barcode = tp.case_barcode AND ss.surgery_number >= s1.surgery_number AND ss.surgery_number < s2.surgery_number)::integer AS received_tmz_sum_cycles,
    (SELECT bool_or(treatment_radiotherapy) FROM clinical.surgeries ss WHERE ss.case_barcode = tp.case_barcode AND ss.surgery_number >= s1.surgery_number AND ss.surgery_number < s2.surgery_number) AS received_rt,
    (SELECT sum(treatment_radiation_dose_gy) FROM clinical.surgeries ss WHERE ss.case_barcode = tp.case_barcode AND ss.surgery_number >= s1.surgery_number AND ss.surgery_number < s2.surgery_number) AS received_rt_sum_gy,
    (SELECT bool_or(treatment_alkylating_agent) FROM clinical.surgeries ss WHERE ss.case_barcode = tp.case_barcode AND ss.surgery_number >= s1.surgery_number AND ss.surgery_number < s2.surgery_number) AS received_alk,
    (SELECT bool_or(treatment_alkylating_agent or treatment_radiotherapy) FROM clinical.surgeries ss WHERE ss.case_barcode = tp.case_barcode AND ss.surgery_number >= s1.surgery_number AND ss.surgery_number < s2.surgery_number) AS received_treatment,
    (CASE
     WHEN r2.grade_int > r1.grade_int THEN 'Grade up'
     WHEN r2.grade_int = r1.grade_int THEN 'Grade stable'
     WHEN r2.grade_int < r1.grade_int THEN 'Grade down'
     END) AS grade_change,
    s2.surgical_interval_mo - s1.surgical_interval_mo AS surgical_interval,
    1 AS surgery,
    cas1.case_overall_survival_mo AS OS,
    CASE WHEN cas1.case_vital_status='alive' THEN 0 WHEN cas1.case_vital_status='dead' THEN 1 END AS case_vital_status,
    (CASE
     WHEN mf2.coverage_adj_mut_freq > mf1.coverage_adj_mut_freq AND mf2.coverage_adj_mut_freq > 10 THEN true 
     WHEN mf2.coverage_adj_mut_freq < 10 THEN false
     ELSE NULL
     END) hypermutator_status,
    (CASE
     WHEN mf2.coverage_adj_mut_freq > mf1.coverage_adj_mut_freq AND mf2.coverage_adj_mut_freq > 10 AND (SELECT bool_or(treatment_alkylating_agent) FROM clinical.surgeries ss WHERE ss.case_barcode = tp.case_barcode AND ss.surgery_number >= s1.surgery_number AND ss.surgery_number < s2.surgery_number) IS TRUE THEN true 
     WHEN mf2.coverage_adj_mut_freq < 10 AND (SELECT bool_or(treatment_alkylating_agent) FROM clinical.surgeries ss WHERE ss.case_barcode = tp.case_barcode AND ss.surgery_number >= s1.surgery_number AND ss.surgery_number < s2.surgery_number) IS TRUE THEN false
     ELSE NULL
     END) alk_assoc_hypermutator_status,
     CASE WHEN nd1.rneo < 1 THEN 1 WHEN nd1.rneo >=1 THEN 0 END AS depletion_initial,
     CASE WHEN nd2.rneo < 1 THEN 1 WHEN nd2.rneo >=1 THEN 0 END AS depletion_recurrent
FROM analysis.tumor_pairs tp
LEFT JOIN biospecimen.aliquots a1 ON a1.aliquot_barcode = tp.tumor_barcode_a
LEFT JOIN biospecimen.aliquots a2 ON a2.aliquot_barcode = tp.tumor_barcode_b
LEFT JOIN clinical.surgeries s1 ON s1.sample_barcode = a1.sample_barcode
LEFT JOIN clinical.surgeries s2 ON s2.sample_barcode = a2.sample_barcode
LEFT JOIN roman_to_int r1 ON r1.grade = s1.grade
LEFT JOIN roman_to_int r2 ON r2.grade = s2.grade
LEFT JOIN analysis.mut_freq mf1 ON mf1.aliquot_barcode = tp.tumor_barcode_a
LEFT JOIN analysis.mut_freq mf2 ON mf2.aliquot_barcode = tp.tumor_barcode_b
LEFT JOIN analysis.neoantigen_depletion nd1 ON nd1.aliquot_barcode = tp.tumor_barcode_a
LEFT JOIN analysis.neoantigen_depletion nd2 ON nd2.aliquot_barcode = tp.tumor_barcode_b
LEFT JOIN clinical.subtypes clin ON clin.case_barcode = s1.case_barcode
LEFT JOIN clinical.cases cas1 ON cas1.case_barcode = s1.case_barcode
INNER JOIN analysis.gold_set gs ON gs.case_barcode = tp.case_barcode
WHERE nd1.rneo IS NOT NULL AND nd2.rneo IS NOT NULL AND (nd1.nobs >= 3 AND nd2.nobs >= 3)
"

res <- dbGetQuery(con, q)

#surg_cox_a <- coxph(Surv(surgical_interval,surgery) ~ nd_a + subtype, data = res)
#surg_cox_b <- coxph(Surv(surgical_interval,surgery) ~ nd_b + subtype, data = res)

os_cox_a <- summary(coxph(Surv(os,case_vital_status) ~ depletion_initial + age + subtype, data = res))
os_cox_b <- summary(coxph(Surv(os,case_vital_status) ~ depletion_recurrent + age + subtype, data = res))

surg_cox_a <- coxph(Surv(surgical_interval,surgery) ~ nd_a + subtype + received_treatment, data = res)
surg_cox_b <- coxph(Surv(surgical_interval,surgery) ~ nd_b + subtype + received_treatment, data = res)


cor(res[,"surgical_interval"],res[,"nd_b"],method="s",use="complete")

#Neoantigen load
#----------------
library(odbc)
library(DBI)
library(ggplot2)
library(survival)

con <- DBI::dbConnect(odbc::odbc(), "VerhaakDB")  

q <- "WITH roman_to_int (grade, grade_int)
AS
(
     VALUES ('I', 1),
            ('II', 2),
            ('III', 3),
            ('IV', 4)
)
SELECT
    tp.tumor_pair_barcode,
    tp.case_barcode,
    tp.tumor_barcode_a,
    tp.tumor_barcode_b,
	clin.idh_codel_subtype AS subtype,
	cas1.case_age_diagnosis_years AS age,
    (CASE
     WHEN s1.surgery_location = s2.surgery_location AND (s1.surgery_laterality = s2.surgery_laterality OR (s1.surgery_laterality IS NULL AND s2.surgery_laterality IS NULL)) THEN 'Local'
     WHEN s1.surgery_location <> s2.surgery_location OR s1.surgery_laterality <> s2.surgery_laterality THEN 'Distal'
     END) AS recurrence_location,
    (SELECT bool_or(treatment_tmz) FROM clinical.surgeries ss WHERE ss.case_barcode = tp.case_barcode AND ss.surgery_number >= s1.surgery_number AND ss.surgery_number < s2.surgery_number) AS received_tmz,
    (SELECT sum(treatment_tmz_cycles) FROM clinical.surgeries ss WHERE ss.case_barcode = tp.case_barcode AND ss.surgery_number >= s1.surgery_number AND ss.surgery_number < s2.surgery_number)::integer AS received_tmz_sum_cycles,
    (SELECT bool_or(treatment_radiotherapy) FROM clinical.surgeries ss WHERE ss.case_barcode = tp.case_barcode AND ss.surgery_number >= s1.surgery_number AND ss.surgery_number < s2.surgery_number) AS received_rt,
    (SELECT sum(treatment_radiation_dose_gy) FROM clinical.surgeries ss WHERE ss.case_barcode = tp.case_barcode AND ss.surgery_number >= s1.surgery_number AND ss.surgery_number < s2.surgery_number) AS received_rt_sum_gy,
    (SELECT bool_or(treatment_alkylating_agent) FROM clinical.surgeries ss WHERE ss.case_barcode = tp.case_barcode AND ss.surgery_number >= s1.surgery_number AND ss.surgery_number < s2.surgery_number) AS received_alk,
    (SELECT bool_or(treatment_alkylating_agent or treatment_radiotherapy) FROM clinical.surgeries ss WHERE ss.case_barcode = tp.case_barcode AND ss.surgery_number >= s1.surgery_number AND ss.surgery_number < s2.surgery_number) AS received_treatment,
    (CASE
     WHEN r2.grade_int > r1.grade_int THEN 'Grade up'
     WHEN r2.grade_int = r1.grade_int THEN 'Grade stable'
     WHEN r2.grade_int < r1.grade_int THEN 'Grade down'
     END) AS grade_change,
    s2.surgical_interval_mo - s1.surgical_interval_mo AS surgical_interval,
    1 AS surgery,
    cas1.case_overall_survival_mo AS OS,
    CASE WHEN cas1.case_vital_status='alive' THEN 0 WHEN cas1.case_vital_status='dead' THEN 1 END AS case_vital_status,
    (CASE
     WHEN mf2.coverage_adj_mut_freq > mf1.coverage_adj_mut_freq AND mf2.coverage_adj_mut_freq > 10 THEN true 
     WHEN mf2.coverage_adj_mut_freq < 10 THEN false
     ELSE NULL
     END) hypermutator_status,
    (CASE
     WHEN mf2.coverage_adj_mut_freq > mf1.coverage_adj_mut_freq AND mf2.coverage_adj_mut_freq > 10 AND (SELECT bool_or(treatment_alkylating_agent) FROM clinical.surgeries ss WHERE ss.case_barcode = tp.case_barcode AND ss.surgery_number >= s1.surgery_number AND ss.surgery_number < s2.surgery_number) IS TRUE THEN true 
     WHEN mf2.coverage_adj_mut_freq < 10 AND (SELECT bool_or(treatment_alkylating_agent) FROM clinical.surgeries ss WHERE ss.case_barcode = tp.case_barcode AND ss.surgery_number >= s1.surgery_number AND ss.surgery_number < s2.surgery_number) IS TRUE THEN false
     ELSE NULL
     END) alk_assoc_hypermutator_status,
     CAST(COALESCE(nd1.neoag_count,0) AS numeric) AS neoag_load_a,
     CAST(COALESCE(nd2.neoag_count,0) AS numeric) AS neoag_load_b
FROM analysis.tumor_pairs tp
LEFT JOIN biospecimen.aliquots a1 ON a1.aliquot_barcode = tp.tumor_barcode_a
LEFT JOIN biospecimen.aliquots a2 ON a2.aliquot_barcode = tp.tumor_barcode_b
LEFT JOIN clinical.surgeries s1 ON s1.sample_barcode = a1.sample_barcode
LEFT JOIN clinical.surgeries s2 ON s2.sample_barcode = a2.sample_barcode
LEFT JOIN roman_to_int r1 ON r1.grade = s1.grade
LEFT JOIN roman_to_int r2 ON r2.grade = s2.grade
LEFT JOIN analysis.mut_freq mf1 ON mf1.aliquot_barcode = tp.tumor_barcode_a
LEFT JOIN analysis.mut_freq mf2 ON mf2.aliquot_barcode = tp.tumor_barcode_b
LEFT JOIN analysis.neoag_freq nd1 ON nd1.aliquot_barcode = tp.tumor_barcode_a
LEFT JOIN analysis.neoag_freq nd2 ON nd2.aliquot_barcode = tp.tumor_barcode_b
LEFT JOIN clinical.subtypes clin ON clin.case_barcode = s1.case_barcode
LEFT JOIN clinical.cases cas1 ON cas1.case_barcode = s1.case_barcode
INNER JOIN analysis.gold_set gs ON gs.case_barcode = tp.case_barcode
"

res <- dbGetQuery(con, q)

os_cox_a <- summary(coxph(Surv(os,case_vital_status) ~ neoag_load_a + age + subtype, data = res))
os_cox_b <- summary(coxph(Surv(os,case_vital_status) ~ neoag_load_b + age+  subtype, data = res))


#Immunotherapy
#----------------


library(odbc)
library(DBI)
library(ggplot2)
library(survival)

con <- DBI::dbConnect(odbc::odbc(), "VerhaakDB")  

q <- "
WITH roman_to_int (grade, grade_int)
AS
(
     VALUES ('I', 1),
            ('II', 2),
            ('III', 3),
            ('IV', 4)
)
SELECT
    tumor_pair_barcode,
    tp.case_barcode,
    tumor_barcode_a,
    tumor_barcode_b,
	clin.idh_codel_subtype AS subtype,
    (CASE
     WHEN s1.surgery_location = s2.surgery_location AND (s1.surgery_laterality = s2.surgery_laterality OR (s1.surgery_laterality IS NULL AND s2.surgery_laterality IS NULL)) THEN 'Local'
     WHEN s1.surgery_location <> s2.surgery_location OR s1.surgery_laterality <> s2.surgery_laterality THEN 'Distal'
     END) AS recurrence_location,
    (SELECT bool_or(treatment_tmz) FROM clinical.surgeries ss WHERE ss.case_barcode = tp.case_barcode AND ss.surgery_number >= s1.surgery_number AND ss.surgery_number < s2.surgery_number) AS received_tmz,
    (SELECT sum(treatment_tmz_cycles) FROM clinical.surgeries ss WHERE ss.case_barcode = tp.case_barcode AND ss.surgery_number >= s1.surgery_number AND ss.surgery_number < s2.surgery_number)::integer AS received_tmz_sum_cycles,
    (SELECT bool_or(treatment_radiotherapy) FROM clinical.surgeries ss WHERE ss.case_barcode = tp.case_barcode AND ss.surgery_number >= s1.surgery_number AND ss.surgery_number < s2.surgery_number) AS received_rt,
    (SELECT sum(treatment_radiation_dose_gy) FROM clinical.surgeries ss WHERE ss.case_barcode = tp.case_barcode AND ss.surgery_number >= s1.surgery_number AND ss.surgery_number < s2.surgery_number) AS received_rt_sum_gy,
    (SELECT bool_or(treatment_alkylating_agent) FROM clinical.surgeries ss WHERE ss.case_barcode = tp.case_barcode AND ss.surgery_number >= s1.surgery_number AND ss.surgery_number < s2.surgery_number) AS received_alk,
    (SELECT bool_or(treatment_alkylating_agent or treatment_radiotherapy) FROM clinical.surgeries ss WHERE ss.case_barcode = tp.case_barcode AND ss.surgery_number >= s1.surgery_number AND ss.surgery_number < s2.surgery_number) AS received_treatment,
    (CASE
     WHEN r2.grade_int > r1.grade_int THEN 'Grade up'
     WHEN r2.grade_int = r1.grade_int THEN 'Grade stable'
     WHEN r2.grade_int < r1.grade_int THEN 'Grade down'
     END) AS grade_change,
    s2.surgical_interval_mo - s1.surgical_interval_mo AS surgical_interval,
    1 AS surgery,
    (CASE
     WHEN mf2.coverage_adj_mut_freq > mf1.coverage_adj_mut_freq AND mf2.coverage_adj_mut_freq > 10 THEN true 
     WHEN mf2.coverage_adj_mut_freq < 10 THEN false
     ELSE NULL
     END) hypermutator_status,
    (CASE
     WHEN mf2.coverage_adj_mut_freq > mf1.coverage_adj_mut_freq AND mf2.coverage_adj_mut_freq > 10 AND (SELECT bool_or(treatment_alkylating_agent) FROM clinical.surgeries ss WHERE ss.case_barcode = tp.case_barcode AND ss.surgery_number >= s1.surgery_number AND ss.surgery_number < s2.surgery_number) IS TRUE THEN true 
     WHEN mf2.coverage_adj_mut_freq < 10 AND (SELECT bool_or(treatment_alkylating_agent) FROM clinical.surgeries ss WHERE ss.case_barcode = tp.case_barcode AND ss.surgery_number >= s1.surgery_number AND ss.surgery_number < s2.surgery_number) IS TRUE THEN false
     ELSE NULL
     END) alk_assoc_hypermutator_status,
	 s1.treatment_chemotherapy_other AS immunotherapy_a,
	 s2.treatment_chemotherapy_other AS immunotherapy_b,	 
     nd1.rneo AS nd_a,
     nd2.rneo AS nd_b
FROM analysis.tumor_pairs tp
LEFT JOIN biospecimen.aliquots a1 ON a1.aliquot_barcode = tp.tumor_barcode_a
LEFT JOIN biospecimen.aliquots a2 ON a2.aliquot_barcode = tp.tumor_barcode_b
LEFT JOIN clinical.surgeries s1 ON s1.sample_barcode = a1.sample_barcode
LEFT JOIN clinical.surgeries s2 ON s2.sample_barcode = a2.sample_barcode
LEFT JOIN roman_to_int r1 ON r1.grade = s1.grade
LEFT JOIN roman_to_int r2 ON r2.grade = s2.grade
LEFT JOIN analysis.mut_freq mf1 ON mf1.aliquot_barcode = tp.tumor_barcode_a
LEFT JOIN analysis.mut_freq mf2 ON mf2.aliquot_barcode = tp.tumor_barcode_b
LEFT JOIN analysis.neoantigen_depletion nd1 ON nd1.aliquot_barcode = tp.tumor_barcode_a
LEFT JOIN analysis.neoantigen_depletion nd2 ON nd2.aliquot_barcode = tp.tumor_barcode_b
LEFT JOIN clinical.subtypes clin ON clin.case_barcode = s1.case_barcode
WHERE (nd1.rneo IS NOT NULL AND nd2.rneo IS NOT NULL AND (nd1.nobs >= 3 AND nd2.nobs >= 3)) AND
	s1.treatment_chemotherapy_other LIKE '%embrolizumab%' OR s2.treatment_chemotherapy_other LIKE '%embrolizumab%' OR
	s1.treatment_chemotherapy_other LIKE 'DC%' OR s2.treatment_chemotherapy_other LIKE 'DC%' OR
	s1.treatment_chemotherapy_other LIKE '%endritic%' OR s2.treatment_chemotherapy_other LIKE '%endritic%' OR
	s1.treatment_chemotherapy_other LIKE '%accine%' OR s2.treatment_chemotherapy_other LIKE '%accine%' OR
	s1.treatment_chemotherapy_other LIKE '%dc%ax%' OR s2.treatment_chemotherapy_other LIKE '%dc%ax%'
"

res <- dbGetQuery(con, q)
