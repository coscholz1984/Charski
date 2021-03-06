# Charski
A colorization tool for point &amp; click pixel-style character art. To run this you need Octave, which is available for many platforms: https://www.gnu.org/software/octave/index , and run the *.m script with it.

![Screenshot of Charski](https://github.com/coscholz1984/Charski/blob/main/Screenshot_v1.jpg?raw=true)

The script has been written and tested in Octave 6.2.0.

A sample for character art with corresponding segmentation maps is provided based on the default template of Adventure Game Studio: https://www.adventuregamestudio.co.uk/ and https://github.com/adventuregamestudio/ags

The tools is based on segmentation maps for individual parts of the image: Hair, skin, eyes, torso, belt, arms, pants, shoes. Individual segments are defined by unique colors (red, green, blue, black, cyan, yellow, black, gray). See the code or example images for exact RGB values. It is possible to colorize individual parts with colors from an extended HSV colormap with two blending modes 'multiply' and 'overlay' (see code for definition, which is equal to the GIMP definition of these blending modes). There is also a gray blending option, which is a kind of multiplicative gamma-filter, but it cannot be combined with color yet. 

Using the delete rows and columns options the sprites can be shrunk in a funny way which gives some more room to modification (specify comma separated integers, corresponding to the coordinates displayed when hovering over the input image). There is also a "growth" option for the hair, which dilates the hair segment of the input frames. Pretty funny, not sure if useful though ;)
