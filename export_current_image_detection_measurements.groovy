import qupath.lib.objects.PathDetectionObject
import qupath.lib.gui.tools.MeasurementExporter
import qupath.lib.objects.PathCellObject


def projectDir = buildFilePath(PROJECT_BASE_DIR)
def exportFolder = buildFilePath(projectDir, "results")
println("results will be saved in :" + exportFolder)
mkdirs(exportFolder)

def imageName = getProjectEntry().getImageName()

// TODO: add specific columns to include instead of all columns
// ex : def columnsToInclude = new String[]{"Area", "Perimeter", "Circularity", "Eccentricity", "Median Intensity", "Nucleus: Area", "Nucleus: Perimeter", ...}

def columnsToInclude = new String[]{""}
def exportPath = buildFilePath(exportFolder, imageName + "_QuPath_detection_measurements.tsv")
def exportFile = new File(exportPath)

// Explicitely target cells instead of all detections

//def exportType = PathDetectionObject.class
def exportType = PathCellObject.class

// println(getProject().getImageList())
print "Results will be saved to project :" + exportPath


def exporter = new MeasurementExporter()
                    .imageList(Arrays.asList(getProjectEntry()))
                    .separator("\t")
                    //.includeOnlyColumns(columnsToInclude) // <-- uncomment for specific list
                    .exportType(exportType)
                    .exportMeasurements(exportFile)
                    
print "Detection measurements export for " + imageName + " complete!"