package molmed.qscripts

import java.io.File

import molmed.utils.{GATKConfig, UppmaxConfig, VariantCallingTarget, VariantCallingUtils}
import org.broadinstitute.gatk.queue.QScript

/**
  * Created by vagrant on 2016-09-21.
  */

case class GenotypeConcordanceTarget(outputDir: File,
                                     compVcf: File,
                                     evalVcf: File
                                     ) extends VariantCallingTarget(outputDir, "baseName", null, Seq.empty[File], None, false, false, 1, Some(compVcf)) {
  override val rawSnpVCF: File = evalVcf
  override val genotypeConcordance: File = new File(outputDir.getPath + "/" + evalVcf.getName.split(".vcf*").head + ".gt_concordance")
}

class GenotypeConcordance  extends QScript {

  /**
    * Command line arguments
    */
  @Input(doc = "one or more files with genotyping results (in vcf format) - the comp set", fullName = "genotypes", shortName = "gt", required = true)
  var snpGenotypes: Seq[File] = Nil

  @Input(doc = "one or more files with variant calls to be checked against the genotypes - the eval set", fullName = "variantcalls", shortName = "vc", required = true)
  var variantCalls: Seq[File] = Nil

  @Argument(doc = "directory where all output will be written", fullName = "outputdir", shortName = "out", required = false)
  var outputDir: String = "genotype_concordance"

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

    val gatkOptions =
      GATKConfig(reference, nbrOfThreads, 1, None, None, None, gatkKey = Option(gatkKey))
    val variantCallingUtils = new VariantCallingUtils(gatkOptions, projectName = Some(UppmaxConfig.defaultProjectId), UppmaxConfig())

    // combine the comp files into one vcf
    val compVcf = new File(outputDir + "/comp_set_combined.vcf.gz")
    this.add(variantCallingUtils.CombineVariantFiles(snpGenotypes, compVcf))

    // do a concordance check for each of the input eval files
    for (target <- variantCalls.map(GenotypeConcordanceTarget(new File(outputDir), compVcf, _))) yield {
      this.add(new variantCallingUtils.SNPGenotypeConcordance(target))
      target.genotypeConcordance
    }

  }

}