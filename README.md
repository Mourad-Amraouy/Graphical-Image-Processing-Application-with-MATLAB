# Graphical-Image-Processing-Application-with-MATLAB
Developed a user-friendly graphical interface application using MATLAB for advanced image processing, enabling efficient manipulation, analysis, and enhancement of digital images.


## Canny filter
The Canny filter is certainly the most known and used filter for edge detection. I will explain step by step the canny filter for contour detection. Step by step because the canny filter is a multi-stage filter. The Canny filter is rarely integrated into a Deep Learning model. So I will describe the different parts and at the same time implement it with Pytorch. It can be customized almost without limit, I have allowed myself some deviations.
![image](https://github.com/Mourad-Amraouy/Graphical-Image-Processing-Application-with-MATLAB/assets/146946535/a6d16852-574b-4567-bf38-12a466084b87)

![image](https://github.com/Mourad-Amraouy/Graphical-Image-Processing-Application-with-MATLAB/assets/146946535/9e51d532-4be8-4f74-ab9a-8e9322bad9de)

## Gaussian filtering
First of all, we usually remove the noise that is present in the input image by applying a blurring filter. The choice of this filter is up to you, but we mostly use a Gaussian filter.
https://miro.medium.com/v2/resize:fit:1100/format:webp/1*13XkXfrnOFC9Mur6Mq2h1A.png

![image](https://github.com/Mourad-Amraouy/Graphical-Image-Processing-Application-with-MATLAB/assets/146946535/c2d075a7-40fa-4c50-850e-706ed993b019)

## Sobel filtering
To detect edges, a filter must be applied to the image to extract the gradients.
                                                                  Sobel kernel on X
                                                                  
![image](https://github.com/Mourad-Amraouy/Graphical-Image-Processing-Application-with-MATLAB/assets/146946535/37d39a4e-5c58-4a84-863c-a6f0898ba2f3)

![image](https://github.com/Mourad-Amraouy/Graphical-Image-Processing-Application-with-MATLAB/assets/146946535/64bbfda5-ea18-40b3-9b35-6cdf3480e1f1)

                                                                  Sobel kernel on Y
The second kernel is used to extract the gradients vertically. The one is the transpose of the other. The two cores have the same role but on different axes.
![image](https://github.com/Mourad-Amraouy/Graphical-Image-Processing-Application-with-MATLAB/assets/146946535/0ec319d6-a213-41f9-914f-9309d8c882d9)
## averaging filter
The averaging filter is an image processing operation used to reduce noise in an image and/or blur an image. For example, applying an averaging filter on the left image gives the right image :

![image](https://github.com/Mourad-Amraouy/Graphical-Image-Processing-Application-with-MATLAB/assets/146946535/71146133-4354-4af1-8e7a-9ca7dd4e15d4) ![image](https://github.com/Mourad-Amraouy/Graphical-Image-Processing-Application-with-MATLAB/assets/146946535/ebd34b08-9732-444b-89ce-35a428c572e9)

The averaging filter belongs to the category of local image filters because to calculate the new value of a pixel, it looks at the value of nearby pixels. Concretely, the filtered value of a pixel p is equal to the average of the values of the pixels close to p.
In general, the "pixels close to p" are defined "as the set of pixels contained in a square of width k centered on p :

![image](https://github.com/Mourad-Amraouy/Graphical-Image-Processing-Application-with-MATLAB/assets/146946535/d43b7de4-79d4-4562-9cf6-370107aeb725)

![image](https://github.com/Mourad-Amraouy/Graphical-Image-Processing-Application-with-MATLAB/assets/146946535/88a2be42-d28a-4f29-92b2-64525518a491)

## Convolution
The convolution, or convolution product, is a generalization of the averaging filter where we consider this time a weighted average. The sliding window is then itself an image which contains the weighting coefficients. It is generally called convolution kernel or convolution mask (kernel or mask in English) :

![image](https://github.com/Mourad-Amraouy/Graphical-Image-Processing-Application-with-MATLAB/assets/146946535/d5dd4fcc-b6c7-4e32-a849-722cddfb0030)

The convolution kernel (in the center) contains the weighting coefficients. The principle is then similar to the averaging filter: to calculate the new value of a pixel on the right, the average of the pixels of the original image (on the left) located under the convolution mask weighted by the mask values is calculated.

## Contour detection
The detection of contours consists in looking for continuous curves along the zones of strong variations in the image. Experiments in neuroscience have indeed shown that the detection of contours is one of the first steps carried out by the visual cortex, suggesting their importance for image analysis processes.
The detection of zones of variation of the gray levels of the image corresponds to the derivation operation. Since a digital image is not a continuous function, the notion of derivative is not formally defined and we will use an analog called gradient. As an image has 2 dimensions, the gradient of the image f, denoted by ∇f, is a vector image, given by the two partial derivatives :
∇f=(∂f∂x,∂f∂y).
In the example below, we can observe an image and its 2 partial derivatives. We note that the horizontal (respectively vertical) variations appear in the partial derivative ∂f∂x(respectively ∂f∂y).

![image](https://github.com/Mourad-Amraouy/Graphical-Image-Processing-Application-with-MATLAB/assets/146946535/21f760ac-6507-4629-9faf-dfd03fdb6712)



## Compute the gradients
Well, we have our gradients on both axes of our image. To detect the contours, we want to have the magnitude of our gradient. We can use the absolute-value norm or the euclidean norm.

![image](https://github.com/Mourad-Amraouy/Graphical-Image-Processing-Application-with-MATLAB/assets/146946535/771cca4a-20ce-4a09-b80b-d1ff3d8e4b18)

![image](https://github.com/Mourad-Amraouy/Graphical-Image-Processing-Application-with-MATLAB/assets/146946535/094b7f2e-4b64-4aa2-b432-78f7d1269ebc)

## Non-Maximum suppression
For thinning the edges, the Non-Maximum Suppression method can be used. Before doing this we need to create the kernels of 45° by 45° directions. (You can refer to this post to understand the rotation matrix)

![image](https://github.com/Mourad-Amraouy/Graphical-Image-Processing-Application-with-MATLAB/assets/146946535/557981ce-d0b3-489d-a46b-b8366aaafa3b)

![image](https://github.com/Mourad-Amraouy/Graphical-Image-Processing-Application-with-MATLAB/assets/146946535/246c4be4-6c0c-4448-afcf-2c6221a19886)

The process will, therefore, require to check the 8-neighborhood (or called Moore’s neighborhood). The concept is fairly simple to understand. For each pixel, we will check the orientation. We are going to see if this pixel is more intense than its neighbor of its gradient’s direction. If yes, then we compare the pixel with its neighbor in the opposite direction. If this pixel has the maximum intensity compared to its two-directional neighbors, then it is the local maximum. This pixel will be kept. In all other cases, it is not a local maximum and the pixel is removed.

![image](https://github.com/Mourad-Amraouy/Graphical-Image-Processing-Application-with-MATLAB/assets/146946535/c214a685-ca32-4513-a782-afa6b092654a)

![image](https://github.com/Mourad-Amraouy/Graphical-Image-Processing-Application-with-MATLAB/assets/146946535/1f25f81f-11b0-4ade-8700-3aaada4d256c)

## Thresholds and Hysteresis
Finally, there only remains to apply thresholds. There are three ways of doing this :

Low-High threshold: the pixels with an intensity higher than the threshold are set to 1 and the others to 0.
Low-Weak and Weak-High thresholds: we set the pixels with high intensity to 1, the pixels with Low intensity to 0 and between the two thresholds we set them to 0.5. They are considered as Weak.
Low-Weak and Weak-High with hysteresis: same as above. Weak pixels are evaluated with their hysteresis level and reassigned as High or Low.

![image](https://github.com/Mourad-Amraouy/Graphical-Image-Processing-Application-with-MATLAB/assets/146946535/c1064a76-3f70-4a7b-b5c6-eae17685ff41)

![image](https://github.com/Mourad-Amraouy/Graphical-Image-Processing-Application-with-MATLAB/assets/146946535/97520c17-f80f-41b9-9654-555ebb3323f5)

![image](https://github.com/Mourad-Amraouy/Graphical-Image-Processing-Application-with-MATLAB/assets/146946535/2ede3da3-bd22-416c-8f47-c74d3e6ad309)


## Mathematical Morphology
The exercises to be performed are located in the code base that you retrieve by registering on the GitHub classroom link received by email 1. Read the repository readme carefully to understand how to use it. The majority of the requested functions already exist in OpenCV: the goal is not to use the OpenCV functions but to code them yourself! We will therefore only use the basic OpenCV containers and the input/output functions.
  ### Median filter
The idea of the median filter takes up the principle of the sliding window but is not interested in a linear combination of the values of the pixels in the window (unlike the convolution product). Instead, we will sort the pixels located under the window in order of increasing gray level and take the pixel located in the middle: we are talking about the median element.
    
For example, if we consider the value list (1,3,6,9,223):
this list is sorted and has 5 elements; the median of this list is the 3rd element, that is, the 6th. Compared to the average, the median has the advantage of being   robust to extreme values which are not necessarily representative of the distribution of values in the list. In the previous case, the average is 48.4, it is strongly influenced by the value 223

![image](https://github.com/Mourad-Amraouy/Graphical-Image-Processing-Application-with-MATLAB/assets/146946535/b4353348-494f-41f9-b0be-6e1cd36111f6)

The median filter is particularly effective for reducing pepper and salt noise where the pixels are randomly transformed into white or black pixels :

![image](https://github.com/Mourad-Amraouy/Graphical-Image-Processing-Application-with-MATLAB/assets/146946535/b3fa338b-a46e-45e7-bf83-2d4d1f4f587f)

  ### Erosion and dilation
  Erosion and dilation are two opposite operations of mathematical morphology. While the filters seen so far focus on modifying the gray levels of pixels, mathematical 
  morphology seeks to modify the shape of objects. Thus, the role of an erosion is to reduce the size of the objects in an image while the dilation will seek to enlarge 
 them.
  #### Dilatation
The expansion operation uses the notion of structuring element which generalizes the concept of sliding window. The idea is that, to be able to study the shape of objects, it is necessary to be able to consider sliding windows with shapes more complex than a square: for example segments, circles, triangles. In practice, a structuring element will simply be a binary image that contains the form of interest.
      
The expansion of a grayscale image f:Z2→R by a binary image SE⊆Z2(the structuring element) is denoted f⊕SEor δSE(f). After expansion, the new value of a pixel is given by the maximum value of the pixels under the structuring element : ∀x,(f⊕SE)(x)=δSE(f)(x)=maxp∈SEf(x+p)
 Compared to the median filter: the window is no longer square and we take the max elementinstead of the median. In order to clearly visualize the result of this operation, it is simpler to limit ourselves to binary images at first :
 
 ![image](https://github.com/Mourad-Amraouy/Graphical-Image-Processing-Application-with-MATLAB/assets/146946535/89eb6a9a-0793-4adb-8d72-20aec1e7e2c5)
 
Example of dilations on binary images. On the left column, we see the original binary image: the black pixels are considered as part of the image. On the middle column, we see the 3 structuring elements: 1 cross, 1 horizontal segment and 1 vertical segment. The right column shows the result of the expansion of the image by the structuring element: the result is a binary image where the black and gray pixels are part of the image (the gray pixels have been added to the left image by the expansion). It can be seen that the expansion enlarges the objects present in the image and that the direction and the size of this enlargement depend on the shape of the structuring element.
In the case of grayscale images, it can be seen that the dilation tends to enlarge the light areas according to the shape of the structuring element :

![image](https://github.com/Mourad-Amraouy/Graphical-Image-Processing-Application-with-MATLAB/assets/146946535/bb1b1fd8-ed05-43f0-aeda-220026624885)

#### Erosion
Erosion is obtained in a similar way to dilation by replacing max per min. Let be a grayscale image f:Z2→R
 and a binary image SE⊆Z2(the structuring element), the erosion of fby SE is denoted by f⊖SE or ϵse(f) and is defined by :
∀x,(f⊕SE)(x)=ϵse(f)(x)=minp∈SEf(x+p)
As for the dilation, in order to clearly visualize the result of this operation, it is simpler to limit ourselves to binary images at first :

![image](https://github.com/Mourad-Amraouy/Graphical-Image-Processing-Application-with-MATLAB/assets/146946535/9e8c5457-c750-45d1-ac47-44728865bd8e)

Example of erosions on binary images. On the left column, we see the original binary image: the black pixels are considered as part of the image. On the middle column, we see the 3 structuring elements: 1 cross, 1 horizontal segment and 1 vertical segment. The right column shows the result of the expansion of the image by the structuring element: the result is a binary image where the black pixels are part of the image (the gray pixels have been removed from the left image by erosion). It can be seen that erosion reduces the size of the objects present in the image and that the direction and size of this narrowing depend on the shape of the structuring element.

![image](https://github.com/Mourad-Amraouy/Graphical-Image-Processing-Application-with-MATLAB/assets/146946535/194e6a05-8096-4fc8-845b-3c2953050674)
![image](https://github.com/Mourad-Amraouy/Graphical-Image-Processing-Application-with-MATLAB/assets/146946535/c1fd0523-d3dd-46dc-bfa2-254595a52cb1)
![image](https://github.com/Mourad-Amraouy/Graphical-Image-Processing-Application-with-MATLAB/assets/146946535/5abab411-3f41-41d8-8cf6-3cd62f24fcfc)

#### Opening and closing
The composition of the expansion and erosion operations makes it possible to produce new operators. The two simplest are closing and opening. For a structuring element is, the ySE closure is defined as an expansion followed by an erosion: ϕSE=ϵSE∘δSE. This operation closes the holes smaller than the structuring element and leaves the rest almost unchanged :

![image](https://github.com/Mourad-Amraouy/Graphical-Image-Processing-Application-with-MATLAB/assets/146946535/f4d7b22b-d6e2-4dd9-a57c-838d889695bf)

#### Morphological gradient
The morphological operators also make it possible to construct a morphological gradient, the morphological counterpart of the notion of linear gradient seen in the chapter on convolution. For an image f and a structuring element is, we will define an internal gradient, an external gradient and a morphological gradient. The examples presented are based on a structuring element in the form of a cross and on the usual binary image :

![image](https://github.com/Mourad-Amraouy/Graphical-Image-Processing-Application-with-MATLAB/assets/146946535/1d836af7-a14d-4ea4-a5ac-959b34d997c4)

The internal gradient is defined by f−(f⊖SE)the difference between the image and its eroded, that is to say the set of pixels removed by erosion. The following result is obtained :

![image](https://github.com/Mourad-Amraouy/Graphical-Image-Processing-Application-with-MATLAB/assets/146946535/21f477fc-4c57-4247-b94d-8b47065e2033)

The external gradient is defined by (f⊕SE)-f the difference between the dilated and the image, that is to say the set of pixels added by the dilation. The following result is obtained :

![image](https://github.com/Mourad-Amraouy/Graphical-Image-Processing-Application-with-MATLAB/assets/146946535/b1bcbb87-3c3a-42a1-bdb4-1e236a3dc89e)

The morphological gradient is defined by (f⊕SE)−(f⊖SE) combination of the two previous approaches.The following result is obtained :

![image](https://github.com/Mourad-Amraouy/Graphical-Image-Processing-Application-with-MATLAB/assets/146946535/9c00c13b-6dd6-4ace-abb4-ab4b6cd2e70c)
The application on grayscale images is done without problems :
![image](https://github.com/Mourad-Amraouy/Graphical-Image-Processing-Application-with-MATLAB/assets/146946535/c57841d1-8648-4a7a-a6e0-cffd478a6ed1)
![image](https://github.com/Mourad-Amraouy/Graphical-Image-Processing-Application-with-MATLAB/assets/146946535/521e6a31-9b33-4559-9a84-3ed1b5b74842)


