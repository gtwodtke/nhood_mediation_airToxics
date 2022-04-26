/*ID SYSTEM DIRECTORIES*/
sysdir

global syspath "C:\Users\wodtke\ado\personal\" 

/*COPY PACKAGE FOLDERS INTO SYS PERSONAL DIRECTORY*/

/*INSTALL PACKAGES FROM LOCAL DISK*/
net install shp2dta, from("${syspath}shp2dta")
net install spmap, from("${syspath}spmap")
net install geodist, from("${syspath}geodist")
