package molmed.qscripts

import molmed.config.FileAndProgramResourceConfig
import molmed.utils._
import molmed.utils.VariantCallingUtils
import molmed.utils.VariantCallingConfig
import org.broadinstitute.gatk.queue.QScript

/**
  * Created by vagrant on 2016-09-21.
  */

class GenotypeConcordance  extends QScript with FileAndProgramResourceConfig {

  /**
    * A helper class for passing VariantCallingTargets to GenotypeConcordance
    * @param outputDir The output directory to write result files to
    * @param compVcf The "truth" set, e.g. a VCF file with genotype calls
    * @param evalVcf The variant calls to compare against the truth set
    */
  case class GenotypeConcordanceTarget(outputDir: File,
                                       compVcf: File,
                                       evalVcf: File
                                      ) extends VariantCallingTarget(outputDir, "baseName", null, Seq.empty[File], None, false, false, 1, Some(compVcf)) {
    override val rawSnpVCF: File = evalVcf
    override val genotypeConcordance: File = new File(outputDir.getPath + "/" + evalVcf.getName.split(".vcf*").head + ".gt_concordance")
  }

  /**
    * Command line arguments
    */
  @Input(doc = "one or more files with genotyping results (in vcf format) - the comp set", fullName = "genotypes", shortName = "gt", required = true)
  var snpGenotypes: Seq[File] = Nil

  @Input(doc = "one or more files with variant calls to be checked against the genotypes - the eval set", fullName = "vcffile", shortName = "vcf", required = false)
  var variantCalls: Seq[File] = Nil

  @Input(doc = "one or more bam files to call variants for and check against the genotypes - the eval set", fullName = "bamfile", shortName = "bam", required = false)
  var bamfiles: Seq[File] = Nil

  @Argument(doc = "directory where all output will be written", fullName = "outputdir", shortName = "out", required = false)
  var outputDir: File = _

  @Argument(doc = "project id for the analysis", fullName = "projectid", shortName = "pid", required = false)
  var projectId: String = UppmaxConfig.defaultProjectId

  @Input(doc = "genome reference file in fasta format", fullName = "reference", shortName = "ref", required = true)
  var reference: File = _

  @Argument(doc = "Number of threads to use by default", fullName = "number_of_threads", shortName = "nt", required = false)
  var nbrOfThreads: Int = 1

  @Argument(doc = "File containing license key for disabling GATK phone home feature", fullName = "gatk_key", shortName = "gatkKey", required = false)
  var gatkKey: File = _

  /**
    * **************************************************************************
    * Main script
    * **************************************************************************
    */

  def script() {

    // verify that we got something to work with
    if ((variantCalls != Nil && bamfiles != Nil) || (variantCalls == Nil && bamfiles == Nil)) {
      throw new IllegalArgumentException("Either variant VCF file(s) or mapped BAM file(s) must be specified")
    }

    // setup configs and utilities
    val gatkOptions =
      GATKConfig(reference, nbrOfThreads, 1, None, None, None, gatkKey = Option(gatkKey))
    val variantCallingUtils = new VariantCallingUtils(gatkOptions, projectName = Some(projectId), UppmaxConfig())
    val variantCallingConfig = VariantCallingConfig(
      qscript = this,
      bams = bamfiles,
      outputDir = outputDir.getOrElse(new File("genotype_concordance")),
      bcftoolsPath = bcftoolsPath,
      runSeparatly = false,
      isExome = false,
      isLowPass = false,
      noRecal = true,
      noIndels = true,
      testMode = false,
      skipAnnotation = true,
      skipVcfCompression = false,
      noBAQ = false)

    // combine the comp files into one vcf
    val compVcf = variantCallingUtils.combineVcfFiles(
      variantCallingConfig, snpGenotypes, new File(projectId + "." + String.valueOf(System.currentTimeMillis()/1000) + ".genotypes.combined.vcf.gz"))

    // if eval VCF files were specified as input, do concordance check on those,
    // else call variants on the bam files and check them
    if (variantCalls != Nil) {
      for (target <- variantCalls.map(GenotypeConcordanceTarget(variantCallingConfig.outputDir, compVcf, _))) yield {
        // create an index for the target file
        val targetIndex = variantCallingUtils.indexVcf(variantCallingConfig, target.evalVcf)
        this.add(variantCallingUtils.GenotypeConcordanceWithIndex(target, targetIndex))
        target.genotypeConcordance
      }
    }
    else if (bamfiles != Nil) {
      // create an updated GATKConfig and VariantCallingUtils with the genotype VCF
      val variantCallingUtilsWithGT = new VariantCallingUtils(
        gatkOptions.copy(snpGenotypingVcf = Some(compVcf)), projectName = Some(projectId), UppmaxConfig())
      variantCallingUtilsWithGT.checkGenotypeConcordance(variantCallingConfig)
    }
  }

}