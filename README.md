# FRANK_V1
FRANK - Resume

FRANK works to optimize the provided structure by attempting to find the global minimum of the potential energy surface. This is achieved using a Gradient Descent Genetic Algorithm which allows for the n-dimensional (n = number of atoms - 6) space in a strategic manner. While the program GAMESS can efficiently run energy calculations and optimize geometries the energy optimization simply traverses down the energy well that the provided structure is in, this often results in simply a local minimum not the global one. 

To overcome this problem the provided structure is translated to an array representing its dihedral angles which makes up each structure’s “genetic code” allowing for point mutations and crossing over to occur to introduce variation in the population. This structures are run through GAMESS to gain fitness scores for individuals and on the outer GGA loop remove individuals found within the same energy wells. 

![unnamed-4](https://user-images.githubusercontent.com/23414761/168407169-59eb7701-1690-4014-bd1e-084987bb02b6.png)
<img width="490" alt="Screen Shot 2022-05-13 at 7 18 20 PM" src="https://user-images.githubusercontent.com/23414761/168407199-820b1c36-470f-404f-89e0-368e9c2fdd78.png">
