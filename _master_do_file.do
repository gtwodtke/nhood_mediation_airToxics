capture clear all 
capture set more off 

global prgmpath "C:\Users\wodtke\Desktop\projects\nhood_mediation_toxins\programs\" 

/*DATA PROCESSING*/
do "${prgmpath}01_create_v01_eclsb.do" nostop

do "${prgmpath}02_create_v02_eclsb_nvars.do" nostop

do "${prgmpath}03_create_v01_ncdb.do" nostop

do "${prgmpath}04_create_v01_caces.do" nostop

do "${prgmpath}05_create_v01_rsei.do" nostop

do "${prgmpath}06_create_v03_eclsb_merged.do" nostop

shell "C:\Program Files\R\R-4.1.3\bin\x64\R.exe" CMD BATCH --vanilla --slave --no-timing --no-echo "${prgmpath}07_create_v04_eclsb_mi.R"
shell DEL "${prgmpath}07_create_v04_eclsb_mi.Rout"
/*
/*ANALYSES*/
do "${prgmpath}08_create_figure_2.do" nostop

shell "C:\Program Files\R\R-4.1.3\bin\x64\R.exe" CMD BATCH --vanilla --slave --no-restore --no-timing --no-echo "${prgmpath}09_create_figure_3.R"
shell DEL "${prgmpath}09_create_figure_3.Rout"

shell "C:\Program Files\R\R-4.1.3\bin\x64\R.exe" CMD BATCH --vanilla --slave --no-restore --no-timing --no-echo "${prgmpath}10_create_figure_4.R"
shell DEL "${prgmpath}10_create_figure_4.Rout"

shell "C:\Program Files\R\R-4.1.3\bin\x64\R.exe" CMD BATCH --vanilla --slave --no-timing --no-echo "${prgmpath}11_create_table_1.R"
shell DEL "${prgmpath}11_create_table_1.Rout"

shell "C:\Program Files\R\R-4.1.3\bin\x64\R.exe" CMD BATCH --vanilla --slave --no-restore --no-timing --no-echo "${prgmpath}12_create_figure_5.R"
shell DEL "${prgmpath}12_create_figure_5.Rout"

shell "C:\Program Files\R\R-4.1.3\bin\x64\R.exe" CMD BATCH --vanilla --slave --no-restore --no-timing --no-echo "${prgmpath}13_create_figure_6.R"
shell DEL "${prgmpath}13_create_figure_6.Rout"

shell "C:\Program Files\R\R-4.1.3\bin\x64\R.exe" CMD BATCH --vanilla --slave --no-restore --no-timing --no-echo "${prgmpath}14_create_figure_7.R"
shell DEL "${prgmpath}14_create_figure_7.Rout"

shell "C:\Program Files\R\R-4.1.3\bin\x64\R.exe" CMD BATCH --vanilla --slave --no-restore --no-timing --no-echo "${prgmpath}15_create_figure_8.R"
shell DEL "${prgmpath}15_create_figure_8.Rout"

shell DEL "${prgmpath}Rplots.pdf"

/*SUPPLEMENTAL MATERIALS*/
do "${prgmpath}16_create_misc_stats.do" nostop

do "${prgmpath}17_create_table_S1.do" nostop

do "${prgmpath}18_create_table_S2.do" nostop

shell "C:\Program Files\R\R-4.1.3\bin\x64\R.exe" CMD BATCH --vanilla --slave --no-restore --no-timing --no-echo "${prgmpath}19_create_figure_S1.R"
shell DEL "${prgmpath}19_create_figure_S1.Rout"

shell "C:\Program Files\R\R-4.1.3\bin\x64\R.exe" CMD BATCH --vanilla --slave --no-restore --no-timing --no-echo "${prgmpath}20_create_figure_S2.R"
shell DEL "${prgmpath}20_create_figure_S2.Rout"

shell "C:\Program Files\R\R-4.1.3\bin\x64\R.exe" CMD BATCH --vanilla --slave --no-restore --no-timing --no-echo "${prgmpath}21_create_figure_S3.R"
shell DEL "${prgmpath}21_create_figure_S3.Rout"

shell "C:\Program Files\R\R-4.1.3\bin\x64\R.exe" CMD BATCH --vanilla --slave --no-restore --no-timing --no-echo "${prgmpath}22_create_figure_S4.R"
shell DEL "${prgmpath}22_create_figure_S4.Rout"


shell "C:\Program Files\R\R-4.1.3\bin\x64\R.exe" CMD BATCH --vanilla --slave --no-restore --no-timing --no-echo "${prgmpath}23_create_table_S3.R"
shell DEL "${prgmpath}23_create_table_S3.Rout"

shell "C:\Program Files\R\R-4.1.3\bin\x64\R.exe" CMD BATCH --vanilla --slave --no-restore --no-timing --no-echo "${prgmpath}24_create_table_S4.R"
shell DEL "${prgmpath}24_create_table_S4.Rout"
*/
shell "C:\Program Files\R\R-4.1.3\bin\x64\R.exe" CMD BATCH --vanilla --slave --no-restore --no-timing --no-echo "${prgmpath}25_create_table_S5.R"
shell DEL "${prgmpath}25_create_table_S5.Rout"

shell "C:\Program Files\R\R-4.1.3\bin\x64\R.exe" CMD BATCH --vanilla --slave --no-restore --no-timing --no-echo "${prgmpath}26_create_table_S6.R"
shell DEL "${prgmpath}26_create_table_S6.Rout"

shell "C:\Program Files\R\R-4.1.3\bin\x64\R.exe" CMD BATCH --vanilla --slave --no-restore --no-timing --no-echo "${prgmpath}27_create_table_S7.R"
shell DEL "${prgmpath}27_create_table_S7.Rout"

shell "C:\Program Files\R\R-4.1.3\bin\x64\R.exe" CMD BATCH --vanilla --slave --no-restore --no-timing --no-echo "${prgmpath}28_create_table_S8.R"
shell DEL "${prgmpath}28_create_table_S8.Rout"

shell "C:\Program Files\R\R-4.1.3\bin\x64\R.exe" CMD BATCH --vanilla --slave --no-restore --no-timing --no-echo "${prgmpath}29_create_table_S9.R"
shell DEL "${prgmpath}29_create_table_S9.Rout"

