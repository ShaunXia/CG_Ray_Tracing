#Project1 | CSC562 Advanced CG 


There are five html files:

 - part0: done 
 - part1: colored_box.htm
 - part2: one_bounce.htm 
 - part3: russian.htm 
 - extra 1: mul_light.htm 
 - extra 3: specular.htm

you can direct run those files in browser and see the results. And in order to make it render faster, I set canvas w*h = 256*256. you can change size in htm files.

To ensure the script is running, you can open console panel to see its progress.

Also I have hard coded my json and js files with my github links in htmls, not in js to avoid corss-domain problem.

    <script type="text/javascript" src="https://pages.github.ncsu.edu/sxia4/cg/triangles.json"></script>

There're some other result images when I wrote this project in result folder.

---
### Finished Points:
**Part 0: Properly turned in assignment**

**Part 1: Add a colored Cornell box, and intersect triangles**

![1][1]

**Part 2: Add one bounce diffuse reflection**
100 Samples with *cosTheta

![2][2]

100 Samples without *cosTheta

![3][3]

20 Samples

![4][4]


**Part 3: Add Russian Roulette**
(p_bounce=0.2)

![5][5]

(p_bounce=0.5)

![6][6]


### Extra Points:
**1% Add multiple lights, and weigh their light appropriately**
Add another light in [0.5,0,0.5]

![7][7]

**1% Add the specular term into your BRDF**
a litter brighter than without specular.

![8][8]


  [1]: https://pages.github.ncsu.edu/sxia4/cg/result/colored_box.png
  [2]: https://pages.github.ncsu.edu/sxia4/cg/result/100sample%20b1_cos_only.png
  [3]: https://pages.github.ncsu.edu/sxia4/cg/result/100sample%20b1_only_add_without_spec.png
  [4]: https://pages.github.ncsu.edu/sxia4/cg/result/sp20_b1_cos.png
  [5]: https://pages.github.ncsu.edu/sxia4/cg/result/100Sample_only_add_with_spec_russian0.8.png
  [6]: https://pages.github.ncsu.edu/sxia4/cg/result/100Sample_only_add_with_spec_russian0.5.png
  [7]: https://pages.github.ncsu.edu/sxia4/cg/result/mul_light.png
  [8]: https://pages.github.ncsu.edu/sxia4/cg/result/100Sample_only_add_with_spec.png