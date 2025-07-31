# blender Seams to Sewing pattern

An add-on for Blender that assists with setting up 2D sewing patterns from 3D models, for cloth simulation and real-life sewing.

![](https://blenderartists.org/uploads/default/optimized/4X/3/7/9/379d4a76a9022a7ff338773500784e22500dd8f6_2_690x207.jpeg)\
[![](https://img.youtube.com/vi/EZr__pTxsKk/mqdefault.jpg)\
▶ Quick guide on youtube](https://www.youtube.com/watch?v=EZr__pTxsKk)

# Goatchurchprime's notes on forked code

This is an extremely well-integrated and interesting plugin.  However it is fundamentally flawed by 
its dependence on the Blender UV unwrapping algorithm that is designed to preserve angles 
not areas because it is intended to make texture mapping look good, not to make 
manufacturable parts.

Luckily there is a function in FreeCAD that does surface flattening properly, and has been used in the 
manufacture of hang-glider designs.  https://forum.freecad.org/viewtopic.php?t=18010

I have bodged the functions enough to use it in Blender with minimal finesse on the UI to design a pattern for a plushie-like object, 
and the results are promising.  Development may continue if there is an interest by someone else in using it as well, 
with the emphasis on a real world application where the accuracy of the pattern matters.

link from the addons directory
> $ ~/.config/blender/4.5/scripts/addons]$ ln -s ../../../../../repositories/blender-seams-to-sewing-pattern-FC/ blender-seams-to-sewing-pattern

go into blender-seams-to-sewing-pattern-FC/magic-blender and execute nix develop to enable python libraries flatmesh and numpy


# installation
Download the archive here:\
https://gitlab.com/thomaskole/blender-seams-to-sewing-pattern/-/archive/master/blender-seams-to-sewing-pattern-master.zip

In Blender, go to `Edit > ⚙️ Preferences > Install..` and select the zip file you just downloaded.\
Enable the add-on in the list.

# very brief overview:
![](https://gitlab.com/thomaskole/blender-seams-to-sewing-pattern/-/wikis/uploads/2364f88e60b43cf0cc44309c2e4f15be/triceratops.gif)

`Object > Seams to Sewing Pattern > Seams to Sewing Pattern`\
turns your mesh into a sewing patten based on it's UV layout.

`Object > Seams to Sewing Pattern > Quick Clothsim`\
Applies some basic cloth sim options to your Object

`Object > Seams to Sewing Pattern > Export Sewing Pattern (.svg)`\
Exports your sewing pattern to a .SVG file for printing and sewing in real life.

`Edge > Clean up Knife Cut`\
Clean up selected edges after you used the knife tool on a mesh

# reporting issues
Something wrong? Please let me know.

There's a Blender Artists thread here: https://blenderartists.org/t/1248713 \
You can also add an issue here on GitLab,\
or get in touch with me otherwise: \
http://www.thomaskole.nl

# troubleshooting

**I'm getting a python error**\
That's bad, please let me know the error and how you triggered it.

**After unfolding, there's some weird long triangles /strips flying out**\
Most likely some non-manifold geometry, overlapping vertices, or bad normals.\
Disable the remesh option and see where it goes wrong in your mesh / UV'seams

**My mesh is imploding on itself during clothsim**\
Yeah, clothsim... Try balancing the "pressure" and "sewing force".\
It can help to keyframe the "pressure" to something very high on frame 1 and decrease over time.

# license

GPL V2:\
[GPL-license.txt](./GPL-license.txt)

# sustainability and ethics 🌱

The production of natural textiles often involve overuse of water, fertilizers and perticides.\
Synthetic textiles are a large, often overlooked source of plastic pollution.\
The textile industry is often associated with unethical work conditions.

If you wish to use this add-on for real-life sewing purposes, I ask of you to use only sustainable / second-hand fabrics, and fair labour.

See also:\
[Sustainability-and-Ethics-notice.txt](./Sustainability-and-Ethics-notice.txt)

# known issues

- UV's are messed up after Remesh
