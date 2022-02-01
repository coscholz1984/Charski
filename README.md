# Charski
A colorization tool for point &amp; click pixel-style character art

The script has been written in octave 4.0.

A sample for character art with corresponding segmentation maps is provided based on the default template of Adventure Game Studio: https://www.adventuregamestudio.co.uk/ and https://github.com/adventuregamestudio/ags

The tools is based on segmentation maps for individual parts of the image: Hair, skin, eyes, torso, belt, arms, pants, shoes. Individual segments are defined by unique colors (red, green, blue, black, cyan, yellow, black, gray). See the code or example images for exact RGB values. It is possible to colorize individual parts with colors from an extended HSV colormap with two blending modes 'multiply' and 'overlay' (see code for definition, which is equal to the GIMP definition of these blending modes). There is also a gray blending option, which is a kind of multiplicative gamma-filter, but it cannot be combined with color yet. 
