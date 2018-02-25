# imageSegmentationRepresentation
Various experimental programs in image processing
The programs are very simplistic trials which may not be of much help to a professional, but can help a student who is attempting to understand or try these concepts.

imageSegmentation.m and objectRepresentation.m are programs for segmenting an image based on thresholding and for representing an object based on object representation techniques like chain code and object signature. Since the images being represented are pedestrians, these techniques are obviously not the right techniques, and should instead be represented by Fourier descriptors.

For the below three programs, accuracy of classification is measured via a Dice score index, since each image has its correspnding ground truth also available. 
segmentEstimator.m is a program that attempts to estimate pixel intensities based on surrounding pixels and hence classify an image into foreground and background. The R squared value is also calculated.
nonLinSegmentEstimator.m is the same as segmentEstimator.m, except that a non linear model is used. The implementation needs to be made more accurate though.
segmentEstimatorHypothesis.m uses the z score to determine whether the hypothesis of a pixel being foreground or background could be rejected or accepted. Done for various alpha values, an ROC curve is also plotted. Also plotted are the sensitivity and specificity.
