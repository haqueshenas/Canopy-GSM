# Canopy-GSM
## Version 2.0
Efficient quantification of intricate shading patterns within 3D vegetation canopies can enhance our understanding of canopy functions and conditions. To quantitatively characterize shading patterns in ground-based images of crop canopies, we developed a simple image mining technique known as the Green-Gradient Based Canopy Segmentation Model (GSM). This model relies on relative variations in RGB triplets under different illuminations.

When running the GSM code, vegetation pixels are categorized into a maximum of 255 groups based on their green levels. Subsequently, the mean red and mean blue levels of each group are calculated and plotted against the green levels. The resulting graph (i.e., the GSM graph) can be readily utilized for:

- Identifying and characterizing canopies, as one or two exponential equations;   
- Classifying canopy pixels based on the degree of exposure to sunlight; and   
- Providing valuable quantitative output values, including the coefficients of the exponential equations and up to 510 coordinates of points from the GSM graph. These outputs can be used in data mining analyses for modern phenotyping.   

**For more detailed information, please refer to the published paper:**   

Haghshenas, A., & Emam, Y. (2020). Green-gradient based canopy segmentation: A multipurpose image mining model with potential use in crop phenotyping and canopy studies. Computers and Electronics in Agriculture, 178, 105740. doi:https://doi.org/10.1016/j.compag.2020.105740  




### Input
>
>  (1) RGB images        
>                              
 
### Output
>
>(1) Processed images
>  - Results of segmentation type 1 (ST1, see the published paper) are saved in the "Processed images" folder. In processed images, the non-vegetation parts (background) are shown in purple.
>
>(2) Graphs
>   - A GSM graph is drawn for each image and saved in the "Graphs" folder.
>
>(3) Tables of results
>> (i) Single .csv files
>>   
For each image, a single .CSV file is made containing various types of GSM outputs (find them in the "Results" folder).
>    These tables are consisted of 22 columns as below: 
>>- "Green level": values 1 to 255 as the range of horizontal axis of GSM graph;
 >>- "Mean red": mean red values at each green level (RGB color system). GSM red curve is drawn using the dataset represented in this column;
 >>- "Green": Exactly the same as the data of "Green level" column; which is repeated here for drawing the 1:1 green line;
 >>- "Mean blue": mean blue values at each green level. GSM blue curve is drawn using the dataset represented in this column;
 >>- "Var. red": variance of red values at each green level;
 >>- "Var. blue": variance of blue values at each green level;
 >>- "Number of pixels": number of pixels contributed in calculations of each green level;
 >>- "Red AUC": the area under the red curve at each green level;
 >>- "Blue AUC": the area under the blue curve at each green level;
 >>- "RB ABC": the area between red and blue curves at each green level;
 >>- "Red slope": local slope of the red curve at each green level;
 >>- "Blue slope": local slope of the blue curve at each green level;
 >>- "Red ST3 class": (result of segmentation type 3;) the ST3 class of the red curve at each green level;
  >>- "Blue ST3 class": the ST3 class of the blue curve at each green level;
  >>- "R2 exp Red": R-squared of fitting the exponential equation to the red curve (see the published paper for the formula);
  >>- "RMSE exp Red": root mean square error of fitting the exponential equation to the red curve;
  >>- "exp-a-red": coefficient "a" of the exponential equation fitted to the red curve;
  >>- "exp-b-red": coefficient "b" of the exponential equation fitted to the red curve;
  >>- "R2 exp blue": R-squared of fitting the exponential equation to the blue curve;
  >>- "RMSE exp blue": root mean square error of fitting the exponential equation to the blue curve;
  >>- "exp-a-blue": coefficient "a" of the exponential equation fitted to the blue curve;
  >>- "exp-b-blue": coefficient "b" of the exponential equation fitted to the blue curve;
  
  >> (ii) Total GSM Exp.csv: this table represents a complete copy of coefficients, R^2, and RMSE of GSM exponential equations stored in all of the single .csv files (i.e. in the 15th to 22nd columns). This file may be useful for datamining approaches.
  
  >> (iii) Total GSM ST2.csv: the total outputs of GSM-ST2 of all single .csv files are collected in this table. The values of red and blue GSM curves (ST2 attributes, i.e. the data of 2nd and 4th columns in each single file) are transposed and stored in a single row.
  The first 255 values in each row represent the data of the red curve (i.e. R_G1 to R_G255), and the second 255 values belong to the blue curve (i.e. B_G1 to B_G255).

## How to run?
 
Step 1: Download and unpack the "CanopyGSM" package.
> - Step 2: Copy and paste your own images into the subfolder "Images".
> - Step 3: Go back to the CanopyGSM directory and open the "main.m" file in Matlab. Run the code!
>          
>
>   Note: 
>   - You can change the segmentation method of ST1 in the "ST1.m" file.
>   -Despite other outputs, results of ST3 (segmentation type 3, columns "Red ST3 class" and "Blue ST3 class" in the individual .csv files) may be exclusive and highly dependent on the canopy status e.g. subjecting to environmental stresses, phenology, etc. Therefore, you may change the two variables of St3_P and St3_m in the code (which are set to 50 and 12 by default), or use alternative methods of curve segmentation. Anyway, changing the ST3 variables, or ignoring this type of segmentation have no effect on other results or the data in other columns. Please see the published paper for more information.

>   - Do not open .csv processors during processing.
>   - If the processing duration is more important for you than having graphs or other additional outputs, you can delete the related code lines. 
>   - If there are many images to be processed using version 2.0, it is recommended to process them in batches of tens or hundreds, depending on the power of your system.
>   - In the output tables, "NaN" stands for Not a Number, which indicates undefined or unrepresentable values (in the present code, it may be resulted by dividing a value by zero, e.g. when there is no pixel with a given green level in the vegetation segment of the image, i.e. Num of pixels=0).
>   - The slopes calculated in the present code are consistent with the outputs of Origin Lab (v. 8).

**Troubleshooting**
In some cases, there may be inconsistencies between the format of input images and the formats defined in the code; even if they are the same. So, if an error occurred during processing, check the upper and lower case letters in file names (e.g. ".jpg" vs ".JPG" or ".tif" vs ".TIF") and edit the code based on the format of your images.


