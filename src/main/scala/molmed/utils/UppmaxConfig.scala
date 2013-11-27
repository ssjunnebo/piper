package molmed.utils

case class UppmaxConfig(
    projId: String = UppmaxConfig.defaultProjectId,
    uppmaxQoSFlag: Option[String] = UppmaxConfig.defaultUppmaxQoSFlag
    )

object UppmaxConfig {
  val defaultProjectId = ""
  val defaultUppmaxQoSFlag: Option[String] = None  
}