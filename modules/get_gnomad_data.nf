def get_vcf_path(version, type, build) {
   switch(version) {
       case ~/4\.1/:
           return "release/${version}/vcf/${type}/gnomad.${type}.v${version}.sites"
       case ~/3\.1\.2/:
           if (build == "GRCh38") {
               return "release/${version}/vcf/${type}/gnomad.${type}.v${version}"
           } else {
               return "release/${version}/vcf/${type}/gnomad.${type}.v${version}.GRCh37"
           }
       case ~/2\.1\.1/:
           return "release/${version}/vcf/${type}/gnomad.${type}.r${version}.sites"
       default:
           error "Unsupported GnomAD version: ${version}"
   }
}

process GET_GNOMAD_DATA {
   tag "${chromosome}:${start}-${end}"
   publishDir "${params.outdir}/gnomad", mode: 'copy'
   
   errorStrategy { task.exitStatus in [255] ? 'ignore' : 'terminate' }

   input:
   tuple val(id), val(chromosome), val(start), val(end)
   val(GNOMAD_release_version)
   val(GNOMAD_data_type)
   val(GNOMAD_genome_build)

   output:
   tuple val(id), val(chromosome), val(start), val(end), path("${chromosome}-${start}-${end}.vcf.bgz"), env(HAS_RECORDS), emit: vcf
   path "versions.yml", emit: versions

   script:
   def base_url = "https://storage.googleapis.com/gcp-public-data--gnomad"
   def vcf_path = get_vcf_path(GNOMAD_release_version, GNOMAD_data_type, GNOMAD_genome_build)
   def vcf_url = "${base_url}/${vcf_path}.${chromosome}.vcf.bgz"
   """
   # Check if file exists using tabix
   if tabix -l ${vcf_url} &>/dev/null; then
       tabix \
           -D \
           --print-header \
           ${vcf_url} \
           "${chromosome}:${start}-${end}" \
       | bcftools view \
           --apply-filters PASS \
           --min-af ${params.min_af} \
           --types snps \
           --threads 1 \
       | bcftools annotate \
           --remove "^INFO/AF,INFO/AC,INFO/AN,INFO/nhomalt" \
           --output-type z \
           --output "${chromosome}-${start}-${end}.vcf.bgz"

       if [ \$(bcftools view -H "${chromosome}-${start}-${end}.vcf.bgz" | head -n 1 | wc -l) -gt 0 ]; then
           HAS_RECORDS="true"
       else
           HAS_RECORDS="false"
       fi
   else
       echo "Error: File does not exist or is not accessible: ${vcf_url}"
       exit 1
   fi

   cat <<-END_VERSIONS > versions.yml
   "${task.process}":
       bcftools: \$(bcftools --version | head -n1 | sed 's/^.*bcftools //; s/ .*\$//')
       tabix: \$(tabix --version | head -n1 | sed 's/^.*tabix //; s/ .*\$//')
   END_VERSIONS
   """
}