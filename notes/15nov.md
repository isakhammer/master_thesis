# 15. November
 Membrane domains and lipid rafts is the same I guess

# Single phase membrane tension dynamics
Modelling and simulations of multi-component lipid
membranes and open membranes via diffuse interface
approaches - Xiaoqiang Wang · Qiang Du

![image](https://user-images.githubusercontent.com/43385748/201888168-883442e7-eeff-496a-9499-d23396a00a90.png)

- Line tension 
"The derived Gamma value of energy per unit length, the line tension g, depends on the elastic moduli of the raft and the surrounding membrane"

Line Tension and Interaction Energies of Membrane Rafts Calculated
from Lipid Splay and Tilt
Peter I. Kuzmin,*y Sergey A. Akimov,*y Yuri A. Chizmadzhev,*y Joshua Zimmerberg,* and Fredric S. Cohen 


# Line tension on lipid rafts related to HIV
Line tension at lipid phase boundaries as driving force for HIV fusion peptide-mediated fusion
Sung-Tae Yang, Volker Kiessling & Lukas K. Tamm

Lipids and proteins are organized in cellular membranes in clusters, often called ‘lipid rafts’. Although raft-constituent ordered lipid domains are thought to be energetically unfavourable for membrane fusion, rafts have long been implicated in many biological fusion processes. For the case of HIV gp41-mediated membrane fusion, this apparent contradiction can be resolved by recognizing that the interfaces between ordered and disordered lipid domains are the predominant sites of fusion. HIV entry into T-cells by membrane fusion

# Common ways to measure surface tension:
https://phys.libretexts.org/Courses/University_of_California_Davis/UCD%3A_Biophysics_241_-_Membrane_Biology/02%3A_Membranes_-_Aggregated_Lipids/2.05%3A_Surface_Tension_and_Line_Tension
![image](https://user-images.githubusercontent.com/43385748/201903583-02e3a111-0a7a-41a7-913d-bdae702ceb57.png)

- High membrane tension hinder endocytosis (fission).  Endocytosis is a cellular process in which substances (vesicles) are brought into the cell. 
[image](https://user-images.githubusercontent.com/43385748/201903791-7d6a0f7a-bb7d-4cbf-8acb-dea39d63752b.png)
 - High membrane tension favors exocytosis (fusion of vesicles). 
![image](https://user-images.githubusercontent.com/43385748/201904434-e16ebf2c-438c-4830-a626-1106b230c6f5.png)

# Multicomponent Lipid Membrane surface tensions dynamics 
 - Modelling and simulations of multi-component lipid
membranes and open membranes via diffuse interface
approaches Xiaoqiang Wang · Qiang Du
https://link.springer.com/content/pdf/10.1007/s00285-007-0118-2.pdf
https://onlinelibrary.wiley.com/doi/pdf/10.1002/pamm.200700446

Assuming two lipids configurations. No phase seperation dynamics, Constant Volume

![image](https://user-images.githubusercontent.com/43385748/201911486-d60ab4af-8d51-482b-b6c7-639722c1c9f3.png)

![image](https://user-images.githubusercontent.com/43385748/201911434-26078266-5c9b-49ee-9250-13fe4822bf80.png)

 -  A SURFACE PHASE FIELD MODEL FOR TWO-PHASE BIOLOGICAL MEMBRANES *CHARLES M. ELLIOTT* AND BJÖRN STINNER
 https://www.jstor.org/stable/pdf/41111185.pdf?refreqid=excelsior%3A86dc81735ef1f028640eeda75aafab67&ab_segments=&origin=&acceptTC=1
Changing out line tension with the Ginzburg-Landau energy (i.e. energy functional to allen cahn/cahn hilliard)

Combining it with the  Canham-Helfrich-Evans energy functional for elasticity for single phase

![image](https://user-images.githubusercontent.com/43385748/201918890-147828f7-a162-4842-a6b3-f684c1804775.png)

s.t. the total energy has the form

![image](https://user-images.githubusercontent.com/43385748/201919039-b73cefeb-af31-4dbd-89a6-c4a8c5cf681e.png)


 - NUMERICAL SHAPE OPTIMIZATION OF THE CANHAM-HELFRICH-EVANS BENDING ENERGY (Ngsolve 2021 )
https://arxiv.org/pdf/2107.13794.pdf
![image](https://user-images.githubusercontent.com/43385748/201913830-cc5c5ccd-df9e-4b64-ab21-8e3918893888.png)



# Variations of the classical Canham Helfrich bending elasticity model

ELASTIC BENDING ENERGY: A VARIATIONAL APPROACH- RICCARDO CAPOVILLA
 -  The paper has alot of good riemann geometry!

The free energy on surface S depends on the shape
functions through the mean curvature K, the intrinsic curvature R, 
![image](https://user-images.githubusercontent.com/43385748/202016337-86282136-9ec1-4a48-8ab1-2a784e6226c5.png)


![image](https://user-images.githubusercontent.com/43385748/202015786-ee4bcd8f-cd5d-4d15-a380-f557dfa495fc.png)

Need to rememebr that if there is a 2D surface with no boundary, then does the Gauss-Bonnet Theorem hold here g is the genus of the surface.  https://en.wikipedia.org/wiki/Gauss%E2%80%93Bonnet_theorem

![image](https://user-images.githubusercontent.com/43385748/202015995-5741b161-22e2-44d2-b017-0781fb42af1c.png)

