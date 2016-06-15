1. Poner los siguientes archivos en el mismo directorio:

MEGA_Consortium_v2_15070954_A1_b138_rsids.txt
list_snps_to_change_alleles.txt
new_positions_snps.txt

2. Actualizar el direccion de este directorio en el script convert_MEGA_alleles_sentli.pl 
 
3. Instrucciones para correr el programa
perl convert_MEGA_alleles.pl archivo.ped 

4. Archivos de salida: 
archivo.flip.ped y archivo.flip.map
duplicate snps archivo.txt: ID de los SNPs duplicados
duplicate positions entrenamiento archivo.txt: ID de los SNPs que tienen la misma posicion
to remove.txt: ID de los SNPs duplicados que se excluyen del archivo map original

5. Otros archivos que se pueden reportar: 
multiallelic SNPs.txt tiene la lista de los SNPs que no son bialelicos (202 SNPs).
complementary snps.txt tiene la lista de SNPs cuyos alelos son C/G o A/T (2002 SNPs).

6. Para correr el script: 
perl convert_MEGA_alleles_sentli_test.pl archivo.ped