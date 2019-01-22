SELECT
    ss.tumor_pair_barcode,
    ss.case_barcode,
    ss.tumor_barcode_a,
    ss.tumor_barcode_b,
    st.idh_codel_subtype,
    snv_driver_count,
    cnv_driver_count,
    snv_driver_count_shared,
    cnv_driver_count_shared,
    snv_driver_count_private_a,
    cnv_driver_count_private_a,
    snv_driver_count_private_b,
    cnv_driver_count_private_b,
    snv_driver_shared,
    cnv_driver_shared,
    snv_driver_stability,
    cnv_driver_stability,
    snv_driver_change,
    cnv_driver_change,
    snv_driver_context_shared,
    cnv_driver_context_shared,
    cnv_driver_context_change,
    snv_driver_context_change,
    snv_driver_evolution
FROM analysis.silver_set ss 
LEFT JOIN analysis.driver_status_snv dss ON ss.tumor_pair_barcode = dss.tumor_pair_barcode
LEFT JOIN analysis.driver_status_cnv dsc ON ss.tumor_pair_barcode = dsc.tumor_pair_barcode
LEFT JOIN clinical.subtypes st ON st.case_barcode = ss.case_barcode