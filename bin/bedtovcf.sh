zcat consensus.norm.vcf.gz | awk '{OFS="\t"; \
	if (!/^#/ && (length($4) > 1 || length($5) > 1))\
	{ print $1,$2-sqrt((length($4)-length($5))^2)-1,$2+sqrt((length($4)-length($5))^2)+1,$4"/"$5,"+" } \
	else if (!/^#/) \
	{ print $1,$2-1,$2,$4"/"$5,"+" } \
	}' | less -S