# CVFinalProject
"Magic" frame project implementation for AIT-Budapest's Computer Vision Applications for Digital Cinema.  
## Usage
path-to-file/imgpro input/0000000.jpg output/0000000.jpg -magic <NUMBER_OF_IMAGES> <FIRST_IMAGE_NUMBER>
Example: src/imgpro input/smallHD/0000118.jpg output/0000118.jpg -magic 60 118
## Step 1: Setup
Iterate through images and call the correct functions on them 
NOTE: Coded based on the name format of jpg files provided by professor (e.g. 0000000.jpg, 0000001.jpg, etc).  
## Step 2: Detect and store the "frame" boundary information
magicFeature()
- Detect feature points (Harris Corner Detector) and determine which corners belong to the trackers on the
corners of the frame. 
- Store the coordinates of the frame corners (center of trackers, calculated by taking the average of the 
coordinates of the trackers' feature points).
NOTE: Tracker detector algorithm based on the trackers that we designed; looks for patches of color around the
feature points and uses RGB color thresholds to determine whether or not the points are part of the tracker).

## Step 3: 
Extract and store the frozen image from the first frame.
magicExtractFrozen()

## Step 4: Local search
Local search in the next frame of the .jpg image sequence to find the new position of the four saved points.

## Step 5: Inverse warping
Calculate inverse transformation between next frame and current frame and replace a portion of the image with the saved image.
magicReplaceFrameContent(nextImage)


## About the skeleton code:
In C++ and provided by the professor. Defined methods read and write jpg images and allow us to overwrite pixel values.
Image processing was written by students.
Professor: Gergely Vass


## Methods added for this project:
**R2Image::magicFeature()** --> feature detection based on assignments from throughout the semester. Called on 
			a single image only.
**R2Image::clusters(x,y)** --> given coordinates (x,y) of a feature point, check if surrounding pixels form clusters of color
			similar to the trackers we created. Return true if clusters of red/blue/green.
**R2Image::findShiftedFrame(nextImage, prev_frame)** --> given frame coordinates prev_frame, returns the new coordinates
			of the frame in the nextImage, running a faster local search to locate the feature points
**R2Image::magicReplaceFrameContent(nextImage, shifted_frame)** --> takes in the next image and replaces the pixels 
			using the coordinates of the shifted frame calculated from findShiftedFrame
**R2Image::magicExtractFrozen()** --> Detects and stores the iformation of the first frame