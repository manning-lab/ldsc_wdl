task preProcess {
	File variant_file

	String? sample_size
	String? number_of_cases
	String? number_of_controls
	String out_file_name
	String? min_info_score
	String? min_maf
	String? no_alleles
	String na_cmd = if defined(no_alleles) then "--no-alleles" else ""
	String? merge_alleles
	String ma_cmd = if defined(merge_alleles) then "--merge-alleles" else ""
	String? min_sample_size
	String? variant_id_column
	String? sample_N_column
	String? sample_case_column
	String? sample_control_column
	String? ref_allele_column
	String? alt_allele_column
	String? p_value_column
	String? frequency_column
	String? signed_summary_stats_column
	String? info_column
	String? nstudy_column
	String? a1_increasing
	String a1_cmd = if defined(a1_increasing) then "--a1-inc" else ""
	String? keep_maf
	String km_cmd = if defined(keep_maf) then "--keep-maf" else ""
	Int memory
	Int disk

	command <<<
		python /ldsc/munge_sumstats.py --sumstats ${variant_file} ${"--N " + sample_size} ${"--N-cas " + number_of_cases} ${"--N-con " + number_of_controls} --out ${out_file_name} ${"--info-min " + min_info_score} ${"--maf-min " + min_maf} ${na_cmd} ${ma_cmd} ${"--n-min " + min_sample_size} ${"--snp " + variant_id_column} ${"--N-col " + sample_N_column} ${"--N-cas-col " + sample_case_column} ${"--N-con-col " + sample_control_column} ${"--a1 " + ref_allele_column} ${"--a2 " + alt_allele_column} ${"--p " + p_value_column} ${"--frq " + frequency_column} ${"--signed-sumstats " + signed_summary_stats_column} ${"--info " + info_column} ${"--nstudy " + nstudy_column} ${a1_cmd} ${km_cmd}
	>>>

	runtime {
		docker: "manninglab/ldsc_wdl:0.1"
		disks: "local-disk ${disk} SSD"
		memory: "${memory} GB"
	}

	output {
		File out_file = "${out_file_name}.sumstats.gz"
		File log = "${out_file_name}.log"
	}
}

task partitionHeritability {
	# add annotation files
	File summary_statistics
	String out_file_name

	Array[File] ld_score_files
	String ld_score_pref = basename(ld_score_files[0],"1.l2.ldscore.gz")
	String ld_score_cmd = "--ref-ld-chr ld/${ld_score_pref}"

	Array[File] number_annotation_files
	Array[File] annotation_files

	Array[File] weight_files
	String weight_pref = basename(weight_files[0],"1.l2.ldscore.gz")
	String weight_cmd = "--w-ld-chr weight/${weight_pref}"

	Array[File] freq_files
	String freq_pref = basename(freq_files[0],"1.frq.gz")
	String freq_cmd = "--frqfile-chr freq/${freq_pref}"

	String? overlap_annotations
	String oa_cmd = if defined(overlap_annotations) then "--overlap-annot" else ""

	Int disk
	Int memory


	command <<<
		mkdir ld && mkdir weight && mkdir freq
		mv -t ld/ ${sep = " " ld_score_files}
		mv -t ld/ ${sep = " " number_annotation_files}
		mv -t ld/ ${sep = " " annotation_files}
		mv -t weight/ ${sep = " " weight_files}
		mv -t freq/ ${sep = " " freq_files}
		python /ldsc/ldsc.py --h2 ${summary_statistics} --out ${out_file_name} ${ld_score_cmd} ${weight_cmd} ${freq_cmd} ${oa_cmd}
	>>>

	runtime {
		docker: "manninglab/ldsc_wdl:0.1"
		disks: "local-disk ${disk} SSD"
		memory: "${memory} GB"
	}

	output {
		File out_file = "${out_file_name}.results"
		File log = "${out_file_name}.log"
	}

}

task plotResults {
	Array[File] result_files
	Array[String] column_labels
	Float? pvalue_threshold

	Int disk
	Int memory

	command <<<
		python /ldsc_wdl/plotLDSC.py --results ${sep = "," result_files} --labels ${sep = "," column_labels} --outpref ${sep = "." column_labels} --pthresh ${default = "0.05", pvalue_threshold}
	>>>

	runtime {
		docker: "manninglab/ldsc_wdl:0.1"
		disks: "local-disk ${disk} SSD"
		memory: "${memory} GB"
	}

	output {
		File enrich = select_first(glob("*enrichment.png"))
		File pvals = select_first(glob("*pvalues.png"))
	}
}

workflow ldsc_wf {
	Array[File] these_variant_files
	String? this_sample_size
	String? this_number_of_cases
	String? this_number_of_controls
	Array[String] these_out_file_names
	String? this_min_info_score
	String? this_min_maf
	String? this_no_alleles
	String? this_merge_alleles
	String? this_min_sample_size
	String? this_variant_id_column
	String? this_sample_N_column
	String? this_sample_case_column
	String? this_sample_control_column
	String? this_ref_allele_column
	String? this_alt_allele_column
	String? this_p_value_column
	String? this_frequency_column
	String? this_signed_summary_stats_column
	String? this_info_column
	String? this_nstudy_column
	String? this_a1_increasing
	String? this_keep_maf
	Int preProcess_memory
	Int preProcess_disk

	Array[File] these_ld_score_files
	Array[File] these_number_annotation_files
	Array[File] these_annotation_files
	Array[File] these_weight_files
	Array[File] these_freq_files
	String? this_overlap_annotations
	Int partitionHeritability_memory
	Int partitionHeritability_disk

	Array[Pair[String,File]] these_names_files = zip(these_out_file_names, these_variant_files)

	scatter (this_name_file in these_names_files) {

		call preProcess {
			input: variant_file = this_name_file.right, sample_size = this_sample_size, number_of_cases = this_number_of_cases, number_of_controls = this_number_of_controls, out_file_name = this_name_file.left, min_info_score = this_min_info_score, min_maf = this_min_maf, no_alleles = this_no_alleles, merge_alleles = this_merge_alleles, min_sample_size = this_min_sample_size, variant_id_column = this_variant_id_column, sample_N_column = this_sample_N_column, sample_case_column = this_sample_case_column, sample_control_column = this_sample_control_column, ref_allele_column = this_ref_allele_column, alt_allele_column = this_alt_allele_column, p_value_column = this_p_value_column, frequency_column = this_frequency_column, signed_summary_stats_column = this_signed_summary_stats_column, info_column = this_info_column, nstudy_column = this_nstudy_column, a1_increasing = this_a1_increasing, keep_maf = this_keep_maf, memory = preProcess_memory, disk = preProcess_disk
		}
		
		call partitionHeritability {
			input: summary_statistics = preProcess.out_file, out_file_name = this_name_file.left, ld_score_files = these_ld_score_files, number_annotation_files = these_number_annotation_files, annotation_files = these_annotation_files, weight_files = these_weight_files, freq_files = these_freq_files, overlap_annotations = this_overlap_annotations, memory = partitionHeritability_memory, disk = partitionHeritability_disk
		}

	}

	call plotResults {
		input: result_files = partitionHeritability.out_file, column_labels = these_out_file_names, memory = partitionHeritability_memory, disk = partitionHeritability_disk
	}

	output {
		Array[File] preProcess_log = preProcess.log
		Array[File] ldsc_log = partitionHeritability.log
		Array[File] ldsc_out = partitionHeritability.out_file
		File heatmap = plotResults.out_png
	}

}