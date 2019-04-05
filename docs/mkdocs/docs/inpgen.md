The FLEUR input generator
===============================

For those, who think that the [Fleur inp.xml](xmlio.md) is too complicated, contains too many options or a too complex format, or  those in need for  defaults for their calculation, a inp-file generator is provided.

The `inpgen` executable takes a simplified input file and generates defaults for:

* the 2D lattice type (squ,p-r,c-r,hex,hx3,obl) and the symmetry information 
* (in film calculations) the vacuum distance and d-tilda
* the atom types and the equivalent atoms
* muffin-tin radii, l-cutoffs and radial grid parameters for the atom types
* plane-wave cutoffs (kmax,gmax,gmaxxc).
* many more specialized parameters ...

In general the input generator does not know:

* is your system magnetic? If some elements (Fe,Co,Ni...) are in the unit cell the program sets `jspins=2` and puts magnetic moments. You might like to change jspins and specify different magnetic moments of our atom types.
* how many k-points will you need? For metallic systems it might be more than for semiconductors or insulators. In the latter cases, also the mixing parameters might be chosen larger.
* is the specified energy range for the valence band ok? Normally, it should, but it's better to check, especially if LO's are used. 

You have to modify your [inp.xml](xmlio.md) file accordingly. Depending on your demands, you might want to change something else, e.g. the XC-functional, the switches for relaxation, use LDA+U etc. ...

#Running inpgen

To call the input generator you tyically do

```
inpgen <simple_file
```

**So please note that the program expects its input from the standard-input.**

The `inpgen` executable accepts a few command-line options. In particular you might find usefull

Option|Description
---|--
`-h`|list off all options
`-explicit`|Some input data that is typically not directly provided in the inp.xml file is now generated. This includes a list of k points and a list of symmetry operations.


#Basic input

Your input should contain (in this order):

* (a) A title
* (b) Optionally: input switches (whether your atomic positions are in internal or scaled Cartesian units)
* (c) Lattice information (either a Bravais-matrix or a lattice type and the required constants (see below); in a.u.)
* (d) Atom information (an identifier (maybe the nuclear number) and the positions (in internal or scaled Cartesian coordinates) ).
* (e) Optionally: for spin-spiral calculations or in case of spin-orbit coupling we need the Q-vector or the Spin-quantization axis to determine the symmetry. 
* (f) Optionally: Preset parameters (Atoms/General)

## Title 
Your title cannot begin with an & and should not contain an ! . Apart from that, you can specify any 80-character string you like.

## Input switches 
The namelist input should start with an &input and end with a / . Possible switches are:

switch|description
---|---
 film=[t,f]       |  if .true., assume film calculation (not necessary if dvac is specified)
 cartesian=[t,f]  |  if .true.,  input is given in scaled Cartesian units,
                  |  if .false., it is assumed to be in internal (lattice) units 
 cal_symm=[t,f]   |  if .true.,  calculate space group symmetry,
                  |  if .false., read in space group symmetry info (file 'sym')
 checkinp=[t,f]   |  if .true.,  program reads input and stops
 inistop=[t,f]    |  if .true.,  program  stops after input file generation (not used now)
 symor=[t,f]      |  if .true.,  largest symmorphic subgroup is selected
 oldfleur=[t,f]   |  if .true.,  only 2D symmetry elements (+I,m_z) are accepted

## An example (including the title): 
```
 3 layer Fe film, p(2x2) unit cell, p4mg reconstruction
 
 &input symor=t oldfleur=t / 
```
## Lattice information 
There are two possibilities to input the lattice information: either you specify the Bravais matrix (plus scaling information) or the Bravais lattice and the required information (axis lengths and angles).

**First variant:**

The first 3 lines give the 3 lattice vectors; they are in scaled Cartesian coordinates. Then an overall scaling factor (aa) is given in a single line and independent (x,y,z) scaling is specified by scale(3) in a following line. For film calculations, the vacuum distance dvac is given in one line together with a3.

Example: tetragonal lattice for a film calculation:
```
  1.0  0.0  0.0        ! a1
  0.0  1.0  0.0        ! a2
  0.0  0.0  1.0  0.9   ! a3 and dvac
  4.89                 ! aa (lattice constant)
  1.0  1.0  1.5        ! scale(1),scale(2),scale(3)
```
The overall scale is set by aa and scale(:) as follows: assume that we want the lattice vectors to be given by
```
 a_i = ( a_i(1) xa , a_i(2) xb , a_i(3) xc ) 
```
then choose aa, scale such that: xa = aa * scale(1)), etc. To make it easy to input sqrts, if scale(i)<0, then scale = sqrt(|scale|) Example: hexagonal lattice
```
      a1 = ( sqrt(3)/2 a , -1/2 a , 0.      )
      a2 = ( sqrt(3)/2 a ,  1/2 a , 0.      )
      a3 = ( 0.          , 0.     , c=1.62a ) 
```
You could specify the following:
```
       0.5  -0.5  0.0     ! a1
       0.5   0.5  0.0     ! a2
       0.0   0.0  1.0     ! a3
       6.2                ! lattice constant
      -3.0   0.0  1.62    ! scale(2) is 1 by default
```

**Second variant:**

Alternatively, you may specify the lattice name and its parameters in a namelist input, e.g.
```
 &lattice latsys='tP' a=4.89 c=6.9155 /
```
The following arguments are implemented: `latsys`, `a0` (default: 1.0), `a`, `b` (default: `a`), `c` (default: `a`), `alpha` (90 degree), `beta` (90 degree), `gamma` (90 degree). Hereby, `latsys` can be chosen from the following table (intended to work for all entries, up to now not all lattices work). `a0` is the overall scaling factor.


 full name                |   No|  short| other_possible_names|   Description            |    Variants
 -------------------------|-----|-------|---------------------|------------------------|---
 simple-cubic             |    1  |cub  | cP, sc | cubic-P
 face-centered-cubic      |    2  |fcc  | cF, fcc| cubic-F
 body-centered-cubic      |    3  |bcc  | cI, bcc| cubic-I
 hexagonal                |    4  |hcp  | hP, hcp| hexagonal-P                   | (15)
 rhombohedral             |    5  |rho  | hr, r,R| hexagonal-R                   | (16)
 simple-tetragonal        |    6  |tet  | tP, st | tetragonal-P
 body-centered-tetragonal |    7  |bct  | tI, bct| tetragonal-I
 simple-orthorhombic      |    8  |orP  | oP     | orthorhombic-P
 face-centered-orthorhombic|   9  |orF  | oF     | orthorhombic-F
 body-centered-orthorhombic|  10  |orI  | oI     | orthorhombic-I
 base-centered-orthorhombic|  11  |orC  | oC, oS | orthorhombic-C, orthorhombic-S |(17,18)
 simple-monoclinic         |  12  |moP  | mP     | monoclinic-P
 centered-monoclinic       |  13  |moC  | mC     | monoclinic-C                   |(19,20)
 triclinic                 |  14  |tcl  | aP|


full name                |   No|  short| other_possible_names|   Description
-------------------------|-----|-------|---------------------|---------------------------
hexagonal2                 | 15|  hdp   |         | hexagonal-2       (60 degree angle)
rhombohedral2              | 16|  trg   |hR2,r2,R2|hexagonal-R2
base-centered-orthorhombic2 |17|  orA   |oA|       orthorhombic-A    (centering on A)
base-centered-orthorhombic3 |18|  orB   |oB|       orthorhombic-B    (centering on B)
centered-monoclinic2        |19|  moA   |mA|       monoclinic-A      (centering on A)
centered-monoclinic3        |20|  moB   |mB|       monoclinic-B      (centering on B)

You should give the independent lattice parameters `a,b,c` and angles `alpha,beta,gamma` as far as required.

## Atom information 

First you give the number of atoms in a single line. If this number is negative, then we assume that only the representative atoms are given; this requires that the space group symmetry be given as input (see below).

Following are, for each atom in a line, the atomic identification number and the position. The identification number is used later as default for the nuclear charge (Z) of the atom. (When all atoms are specified and the symmetry has to be found, the program will try to relate all atoms of the same identifier by symmetry. If you want to manipulate specific atoms later (e.g. change the spin-quantization axis) you have to give these atoms different identifiers. Since they can be non-integer, you can e.g. specify 26.01 and 26.02 for two inequivalent Fe atoms, only the integer part will be used as Z of the atom.)

The input of the atomic positions can be either in scaled Cartesian or lattice vector units, as determined by logical `cartesian` (see above). For supercells, sometimes more natural to input positions in scaled Cartesian.

A possible input (for  CsCl ) would be:
```
  2
 55 0.0 0.0 0.0
 17 0.5 0.5 0.5 
```
or, for a p4g reconstructed Fe trilayer specifying the symmetry:
```
 -2
 26 0.00 0.00 0.0 
 26 0.18 0.32 2.5 

 &gen         3

   -1    0    0        0.00000
    0   -1    0        0.00000
    0    0   -1        0.00000

    0   -1    0        0.00000
    1    0    0        0.00000
    0    0    1        0.00000

   -1    0    0        0.50000
    0    1    0        0.50000
    0    0    1        0.00000 /
```
Here, `&gen` indicates, that just the generators are listed (the 3×3 block is the rotation matrix [only integers], the floating numbers denote the shift); if you like to specify all symmetry elements, you should start with `&sym`. You have furthermore the possibility to specify a global shift of coordinates with e.g.
```
 &shift 0.5 0.5 0.5 /
```
or, to introduce additional scaling factors
```
 &factor 3.0 3.0 1.0 /
```
by which your atomic positions will be divided (the name "factor" is thus slightly counterintuitive).

## Ending an input file

If inpgen.x should stop reading the file earlier (e.g. you have some comments below in the file) or if inpgen.x fails to recognize the end of the input file (which happens with some compilers), one can use the following line:

 `&end /`

## Special cases 

### Film calculations 

In the case of a film calculation, the surface normal is always chosen in z-direction. A two-dimensional Bravais lattice correspond then to the three-dimensional one according to the following table:

 lattice|description
 --|--
 square                    |  primitive tetragonal 
 primitive rectangular     |  primitive orthorhombic
 centered rectangular      |  base centered orthorhombic
 hexagonal                 |  hexagonal
 oblique                   |  monoclinic

The z-coordinates of all atoms have to be specified in Cartesian units (a.u.), since there is no lattice in the third dimension, to which these values could be referred. Since the vacuum boundaries will be chosen symmetrically around the z=0 plane (i.e. -dvac/2 and dvac/2), the  atoms should also be placed symmetrically around this plane.
   
The initial values specified for `a3` and `dvac` (i.e. the third dimension, see above) will be adjusted automatically so that all atoms fit in the unit cell. This only works if the atoms have been placed symmetrically around the z=0 plane.

### Spin-spiral or SOC 

If you intend a spin-spiral calculation, or to include spin-orbit coupling, this can affect the symmetry of your system:

* a spin spiral in the direction of some vector q is only consistent with symmetry elements that operate on a plane perpendicular to q, while
* (self-consistent) inclusion of spin-orbit coupling is incompatible with any symmetry element that changes the polar vector of orbital momentum L that is supposed to be parallel to the spin-quantization axis (SQA) 

Therefore, we need to specify either q or the SQA, e.g.:
```
 &qss 0.0 0.0 0.1 /
```
(the 3 numbers are the x,y,z components of q) to specify a spin-spiral in z-direction, or
```
 &soc 1.5708 0.0 /
```
(the 2 numbers are theta and phi of the SQA) to indicate that the SQA is in x-direction.

Be careful if symmetry operations that are compatible with the chosen q-vector relate two
atoms in your structure, they also will have the same SQA in the muffin-tins!   

## Preset parameters

### Atoms

After you have given general information on your setup, you can specify a number of parameters for one or several atoms that the input-file generator will use while generating the inp file instead of determining the values by itself. The list of parameters for one atom must contain a leading `&atom` flag and end with a `/`. You have to specify the atom for which you set the parameters by using the parameter `element`. If there are more atoms of the same element, you can specify the atom you wish to modify by additionally setting the `id` tag. All parameters available are

parameter|description
--|--
 id=[atomic identification number]      |   identifies the atom you wish to modify.
 z=[charge number]                      |   specifies the charge number of the atom.
 rmt=[muffin-tin radius]                |   specifies a muffin-tin radius for the atom to modify.
 dx=[log increment]                     |   specifies the logarithmic increment of the radial mesh for the atom to modify.
 jri=[# mesh points]                    |   specifies the number of mesh points of the radial mesh for the atom to modify.
 lmax=[spherical angular momentum]      |   specifies the maximal spherical angular momentum of the atom to modify.
 lnonsph=[nonspherical angular momentum]|   specifies the maximal angular momentum up to which non-spherical parts are included to quantities of the atom to modify.
 ncst=[number of core state]            |   specifies the number of states you wish to include in the core of the atom to modify.
 econfig=[core states\|valence states]   |   specifies, which states of the atom to modify are put into the core and which are treated as valence states. This is a string.  You can use `[element name of noble gas]` to shorten the list. d and f states will be filled preferring magnetization.
 bmu=[magnetic moment]                   |  specifies the magnetic moment of the atom to modify.
 lo=[list of local orbitals]             |  specifies, which states shall be treated as local orbitals. This is a string.
 element=[name of the element]           |  identifies the atom to modify by its element name. This is a string. You must specify this.

### General

You also might want to set more general parameters like the choice of the exchange-correlation potential or the desired reciprocal grid in the Brillouin zone beforehand. Those parameters can be given as a namelist using the `&comp`, `&exco`, `&film` and/or `&kpt` flag. The corresponding line in the input-file for the input-file generator has to end with a `/`. All parameters available are, sorted by their affiliation

**&comp**:

parameter|description
--|--
 jspins=[number of spins]              |    specifies the number of spins for the calculation.
 frcor=[frozen core?]                  |    specifies whether or not the frozen-core approximation is used.
 ctail=[core-tail correction?]         |    specifies whether or not the core-tail correction is used.
 kcrel=[fully-magnetic dirac core?]    |    specifies whether or not the core is treated fully-relativistic.
 gmax=[dop PW-cutoff]                  |    specifies the plane-wave cutoff for the density and potential Fourier expansion.
 gmaxxc=[xc-pot PW-cutoff]             |    specifies the plane-wave cutoff for the exchange-correlation potential Fourier expansion.
 kmax=[basis set size]                 |    specifies the cutoff up to which plane-waves are included into the basis.


**&exco**:

parameter|description
--|-
xctyp=[xc-potential]                    |  specifies the choice of the exchange-correlation potential. This is a string.
 relxc=[relativistic?]                  |   specifies whether or not relativistic corrections are used.

**&film**:

parameter|description
--|--
 dvac=[vacuum boundary]                  |  specifies the vacuum boundary in case of film calculations.
 dtild=[z-boundary for 3D-PW box]        |  specifies the z-boundary for the 3D plane-wave box in case of film calculations.

**&kpt**:

parameter|description
--|--
 nkpt=[number of k-pts]                   | specifies the number of k-points in the IBZ to be used.
 div1=[number of k-pts x-direction]       | specifies the exact number of k-points to be used in the full BZ zone along x-direction (equidistant mesh).
 div2=[number of k-pts y-direction]       | specifies the exact number of k-points to be used in the full BZ zone along y-direction (equidistant mesh).
 div3=[number of k-pts z-direction]       | specifies the exact number of k-points to be used in the full BZ zone along z-direction (equidistant mesh).
 tkb=[smearing parameter]                 | specifies a smearing parameter for Gauss- or Fermi-smearing method.
 tria=[triangular method]                | specifies whether or not triangular method shall be used.


Here is an example of an input file for the input file generator of Europium Titanate in which local orbitals are used for the 5s and 5p states of Europium and the muffin-tin radius of one Oxygen atom is manually set. Also, the exchange-correlation potential is chosen to be that of Vosko, Wilk, and Nusair and a k-point mesh is defined for the Brillouin-zone.


````
Europium Titanate Perovskite Structure

 &input cartesian=t inistop=t oldfleur=f /

 &lattice latsys='sc' a= 7.38 a0= 1.0 /

   5 
   63     0.000  0.000  0.000
   08.01  0.000  0.500  0.500
   08.02  0.500  0.000  0.500
   08.03  0.500  0.500  0.000
   22     0.500  0.500  0.500

 &atom element="eu" lo="5s 5p" econfig="[Kr] 4d10|4f7 5s2 5p6 6s2" /
 &atom element="o" id=08.03 rmt=1.17 /
 &exco xctyp='vwn' /
 &kpt div1=5 div2=5 div3=5 /
````