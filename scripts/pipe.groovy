import qupath.ext.instanseg.core.InstanSeg
import qupath.lib.images.servers.ImageServerMetadata
import qupath.lib.images.servers.ImageChannel
import java.io.File

// ========== CONFIGURATION ==========
def modelPath = '/Volumes/ext_TD/01-project/project-life/downloaded/fluorescence_nuclei_and_cells-0.1.1'
def classifierDir = '/Volumes/ext_TD/01-project/project-life/classifiers/object_classifiers'
def classifierMap = [
    'PanCK': 'PanCK_classifier.json',
    'FOXP3': 'FOXP3_classifier.json',
    'CD8': 'CD8_classifier.json',
    'CD163': 'CD163_classifier.json',
    'CAF': 'CAF_classifier.json'
]

// Define channel names
def channelNames = ['DAPI', 'FOXP3', 'PanCK', 'CD8', 'CD163', 'α-SMA']

// Output directory for CSV files
def outputDir = buildFilePath(PROJECT_BASE_DIR, 'results')
mkdirs(outputDir)

// ========== PROCESS CURRENT IMAGE ==========
def imageData = getCurrentImageData()
def hierarchy = imageData.getHierarchy()
def imageName = imageData.getServer().getMetadata().getName()

println "========================================="
println "Processing: ${imageName}"
println "========================================="

// STEP 0: Rename channels FIRST
println "0. Renaming channels..."
def server = imageData.getServer()
def originalMetadata = server.getMetadata()

def newChannels = []
channelNames.eachWithIndex { name, i ->
    if (i < originalMetadata.getChannels().size()) {
        def oldChannel = originalMetadata.getChannels()[i]
        newChannels.add(ImageChannel.getInstance(name, oldChannel.getColor()))
        println "   Channel ${i + 1}: ${name}"
    }
}

def newMetadata = new ImageServerMetadata.Builder(originalMetadata)
    .channels(newChannels)
    .build()

imageData.updateServerMetadata(newMetadata)
println "   ✓ Channels renamed"

// Set image type
setImageType('FLUORESCENCE')

// Clear existing objects
clearAllObjects()

// Step 1: Create full frame annotation
println "1. Creating full frame annotation..."
createFullImageAnnotation(true)

// Step 2: Run InstaSeg
println "2. Running InstaSeg..."
InstanSeg.builder()
    .modelPath(modelPath)
    .device("cpu")
    .nThreads(2)
    .tileDims(256)
    .interTilePadding(16)
    .makeMeasurements(true)
    .randomColors(false)
    .build()
    .detectObjects(imageData)

def detections = getCellObjects()
println "   ✓ Detected ${detections.size()} cells"

if (detections.isEmpty()) {
    println "   ⚠ WARNING: No cells detected!"
    return
}

// Step 3: Apply classifiers sequentially
println "3. Applying classifiers sequentially..."

def cellResults = [:]
detections.each { cell ->
    cellResults[cell] = [:]
}

classifierMap.each { markerName, classifierFile ->
    println "   → Applying ${markerName}..."
    
    def classifier = loadObjectClassifier(new File(classifierDir, classifierFile).getAbsolutePath())
    classifier.classifyObjects(imageData, detections, true)
    
    detections.each { cell ->
        def classification = cell.getPathClass()?.getName()
        def isPositive = (classification == "${markerName}_Positive")
        cellResults[cell][markerName] = isPositive
    }
    
    def positiveCount = cellResults.values().count { it[markerName] }
    println "     ✓ ${markerName} saved (${positiveCount} positive, ${detections.size() - positiveCount} negative)"
}

// Step 4: Combine all results
println "4. Creating combined classifications..."
detections.each { cell ->
    def positiveMarkers = []
    
    classifierMap.keySet().each { marker ->
        if (cellResults[cell][marker]) {
            positiveMarkers.add(marker)
        }
    }
    
    if (positiveMarkers.isEmpty()) {
        cell.setPathClass(getPathClass("Negative"))
    } else {
        cell.setPathClass(getPathClass(positiveMarkers.join("+") + "+"))
    }
}

// Step 5: Export CSV
println "5. Exporting to CSV..."
def csvPath = buildFilePath(outputDir, imageName + '_results.csv')

new File(csvPath).withWriter { writer ->
    writer.writeLine("filename,x,y,class")
    
    detections.each { detection ->
        def roi = detection.getROI()
        def x = roi.getCentroidX()
        def y = roi.getCentroidY()
        def className = detection.getPathClass()?.getName() ?: "Negative"
        
        writer.writeLine("${imageName},${x},${y},${className}")
    }
}

println "   ✓ Saved to: ${csvPath}"

// Show summary
def classCounts = detections.groupBy { 
    it.getPathClass()?.getName() ?: "Negative" 
}.collectEntries { k, v -> [k, v.size()] }

println "\n========================================="
println "Summary:"
classCounts.each { className, count ->
    println "  ${className}: ${count}"
}
println "========================================="
println "COMPLETE!"