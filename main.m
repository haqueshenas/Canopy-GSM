% Canopy GSM (version 2.0)

%   A code for:
%   Green-gradient based canopy segmentation: A multipurpose image mining model with potential use in crop phenotyping and canopy studies
 
  
%    Authors: Abbas Haghshenas (a), Yahya Emam (b), Saeid Jafarizadeh (C)

%    (a) PhD, Department of Plant Production and Genetics, Shiraz University, Shiraz, Iran. 

%    (b) Professor, Department of Plant Production and Genetics, Shiraz University, Shiraz, Iran.

%    (c) Former MSc student, Vision Lab, Electrical and Computer Engineering School, Shiraz University, Shiraz, Iran.
    
%   Email: 
%   abbas.haghshenas@shirazu.ac.ir
%   Yaemam@shirazu.ac.ir 
           
  
%%  The code is used in the paper entitled:
 
%   Green-gradient based canopy segmentation: A multipurpose image mining model with potential use in crop phenotyping and canopy studies
  
%   (Haghshenas, A., & Emam, Y. (2020). Green-gradient based canopy segmentation: A multipurpose image mining model with potential use in crop
%    phenotyping and canopy studies. Computers and Electronics in Agriculture, 178, 105740. doi:https://doi.org/10.1016/j.compag.2020.105740)
%---------------------------------------------------------------------------------------------------------------


% Canopy GSM (version 2.0)

clc
clear
close all

% Setting ST3 variables

% (ST3_P and ST3_m may be set to custom values, depended on the purposes
% and conditions; e.g. in our experiment, they were set to 50 and 20, respectively.)

St3_P = 50;
St3_m = 12;

% Making output subfolders

newSubFolder = sprintf('Graphs', '.\code');
if ~exist(newSubFolder, 'dir')
  mkdir(newSubFolder);
end

newSubFolder = sprintf('Processed images', '.\code');
if ~exist(newSubFolder, 'dir')
  mkdir(newSubFolder);
end

newSubFolder = sprintf('Results', '.\code');
if ~exist(newSubFolder, 'dir')
  mkdir(newSubFolder);
end

% Loading images

disp(['Loading Images from directory (Images)... ']); 
imgpath = './Images/';
imgtype='*.jpg';
images = dir([imgpath imgtype]);
imgtype='*.jpeg';
images = [images ; dir([imgpath imgtype]);]
imgtype='*.tiff';
images = [images ; dir([imgpath imgtype]);]
imgtype='*.tif';
images = [images ; dir([imgpath imgtype]);]

Data_total_st2 = {'Image','R G1','R G2','R G3','R G4','R G5','R G6','R G7','R G8','R G9','R G10','R G11','R G12','R G13',...
'R G14','R G15','R G16','R G17','R G18','R G19','R G20','R G21','R G22','R G23','R G24','R G25',...
'R G26','R G27','R G28','R G29','R G30','R G31','R G32','R G33','R G34','R G35','R G36','R G37',...
'R G38','R G39','R G40','R G41','R G42','R G43','R G44','R G45','R G46','R G47','R G48','R G49',...
'R G50','R G51','R G52','R G53','R G54','R G55','R G56','R G57','R G58','R G59','R G60','R G61',...
'R G62','R G63','R G64','R G65','R G66','R G67','R G68','R G69','R G70','R G71','R G72','R G73',...
'R G74','R G75','R G76','R G77','R G78','R G79','R G80','R G81','R G82','R G83','R G84','R G85',...
'R G86','R G87','R G88','R G89','R G90','R G91','R G92','R G93','R G94','R G95','R G96','R G97',...
'R G98','R G99','R G100','R G101','R G102','R G103','R G104','R G105','R G106','R G107','R G108',...
'R G109','R G110','R G111','R G112','R G113','R G114','R G115','R G116','R G117','R G118','R G119',...
'R G120','R G121','R G122','R G123','R G124','R G125','R G126','R G127','R G128','R G129','R G130',...
'R G131','R G132','R G133','R G134','R G135','R G136','R G137','R G138','R G139','R G140','R G141',...
'R G142','R G143','R G144','R G145','R G146','R G147','R G148','R G149','R G150','R G151','R G152',...
'R G153','R G154','R G155','R G156','R G157','R G158','R G159','R G160','R G161','R G162','R G163',...
'R G164','R G165','R G166','R G167','R G168','R G169','R G170','R G171','R G172','R G173','R G174',...
'R G175','R G176','R G177','R G178','R G179','R G180','R G181','R G182','R G183','R G184','R G185',...
'R G186','R G187','R G188','R G189','R G190','R G191','R G192','R G193','R G194','R G195','R G196',...
'R G197','R G198','R G199','R G200','R G201','R G202','R G203','R G204','R G205','R G206','R G207',...
'R G208','R G209','R G210','R G211','R G212','R G213','R G214','R G215','R G216','R G217','R G218',...
'R G219','R G220','R G221','R G222','R G223','R G224','R G225','R G226','R G227','R G228','R G229',...
'R G230','R G231','R G232','R G233','R G234','R G235','R G236','R G237','R G238','R G239','R G240',...
'R G241','R G242','R G243','R G244','R G245','R G246','R G247','R G248','R G249','R G250','R G251',...
'R G252','R G253','R G254','R G255','B G1','B G2','B G3','B G4','B G5','B G6','B G7','B G8','B G9',...
'B G10','B G11','B G12','B G13','B G14','B G15','B G16','B G17','B G18','B G19','B G20','B G21',...
'B G22','B G23','B G24','B G25','B G26','B G27','B G28','B G29','B G30','B G31','B G32','B G33',...
'B G34','B G35','B G36','B G37','B G38','B G39','B G40','B G41','B G42','B G43','B G44','B G45',...
'B G46','B G47','B G48','B G49','B G50','B G51','B G52','B G53','B G54','B G55','B G56','B G57',...
'B G58','B G59','B G60','B G61','B G62','B G63','B G64','B G65','B G66','B G67','B G68','B G69',...
'B G70','B G71','B G72','B G73','B G74','B G75','B G76','B G77','B G78','B G79','B G80','B G81',...
'B G82','B G83','B G84','B G85','B G86','B G87','B G88','B G89','B G90','B G91','B G92','B G93',...
'B G94','B G95','B G96','B G97','B G98','B G99','B G100','B G101','B G102','B G103','B G104',...
'B G105','B G106','B G107','B G108','B G109','B G110','B G111','B G112','B G113','B G114','B G115',...
'B G116','B G117','B G118','B G119','B G120','B G121','B G122','B G123','B G124','B G125','B G126',...
'B G127','B G128','B G129','B G130','B G131','B G132','B G133','B G134','B G135','B G136','B G137',...
'B G138','B G139','B G140','B G141','B G142','B G143','B G144','B G145','B G146','B G147','B G148',...
'B G149','B G150','B G151','B G152','B G153','B G154','B G155','B G156','B G157','B G158','B G159',...
'B G160','B G161','B G162','B G163','B G164','B G165','B G166','B G167','B G168','B G169','B G170',...
'B G171','B G172','B G173','B G174','B G175','B G176','B G177','B G178','B G179','B G180','B G181',...
'B G182','B G183','B G184','B G185','B G186','B G187','B G188','B G189','B G190','B G191','B G192',...
'B G193','B G194','B G195','B G196','B G197','B G198','B G199','B G200','B G201','B G202','B G203',...
'B G204','B G205','B G206','B G207','B G208','B G209','B G210','B G211','B G212','B G213','B G214',...
'B G215','B G216','B G217','B G218','B G219','B G220','B G221','B G222','B G223','B G224','B G225',...
'B G226','B G227','B G228','B G229','B G230','B G231','B G232','B G233','B G234','B G235','B G236',...
'B G237','B G238','B G239','B G240','B G241','B G242','B G243','B G244','B G245','B G246','B G247',...
'B G248','B G249','B G250','B G251','B G252','B G253','B G254','B G255'};
Data_total_Exp = {'Image','R2_exp_Red','RMSE_exp_Red','exp_a_Red','exp_b_Red','R2_exp_Blue','RMSE_exp_Blue','Exp_a_Blue','Exp_b_Blue'};



for i=2:length(images)+1
Data = {'Green level','Mean Red','Green','Mean Blue',...
'Var Red','Var Blue','Number of pixels','Red_AUC','Blue_AUC','RB_ABC',...
'Red Slope','Blue slope','Red ST3 class','Blue ST3 class','R^2 exp Red','RMSE exp Red ','exp a Red','exp b Red',...
    'R^2 exp Blue','RMSE exp Blue ','exp a Blue','exp b Blue'};
Data4 = {'Image name'};
Image = imread(['./Images/',images(i-1).name]);
R =  (Image(:,:,1));
G =  (Image(:,:,2));
B =  (Image(:,:,3));
disp('=================== New image =======================')
disp(['Processing image(',num2str(i-1),':',images(i-1).name,')... ']);

ImageForgraoundRGB = Image;

BinaryImageRGB = ST1(R,G,B);


ImageForgraound =[];
ImageForgraound(:,:,1) =uint8(ST1(R,G,B)).*Image(:,:,1);
ImageForgraound(:,:,2) =uint8(ST1(R,G,B)).*Image(:,:,2);
ImageForgraound(:,:,3) =uint8(ST1(R,G,B)).*Image(:,:,3);

% Setting background color

r = 150;
g = 0;
b = 150;
ImageForgraound2 = uint8(ChangeBackColor(ImageForgraound,r,g,b));
imwrite(ImageForgraound2,['./Processed images/Processed_',[images(i-1).name],'.jpg'])
 
R1 = R.*uint8(ST1(R,G,B));
G1 = G.*uint8(ST1(R,G,B));
B1 = B.*uint8(ST1(R,G,B));

% GSM (ST2)

for j=1:255
    disp(['Computing level ',num2str(j),' of ',images(i-1).name,'... ']);

Data(j+2,1) = num2cell(j);
Data(j+2,2) = num2cell(mean(mean(R(find(G1==j)))));
Data(j+2,3) = num2cell(j);
Data(j+2,4) = num2cell(mean(mean(B(find(G1==j)))));

%% variance 
R = double(R);
G = double(G);
B = double(B);

Data(j+2,5) = num2cell((var(R(find(G1==j )))));
Data(j+2,6) = num2cell((var(B(find(G1==j )))));
Data(j+2,7) = num2cell(length(R(find(G1==j ))));

tempp = cell2mat(Data(3:end,2));
tempp(isnan(tempp))= 0;
tempp = sum(tempp);

Data(j+2,8) = num2cell(tempp);
tempp = cell2mat(Data(3:end,4));
tempp(isnan(tempp))= 0;
tempp = sum(tempp);

Data(j+2,9) = num2cell(tempp);
Data(j+2,10) = num2cell(cell2mat(Data(j+2,8))-cell2mat(Data(j+2,9)));

end

for j=1:257

    % Curve derivation
    
    if(~isnan(cell2mat(Data(j,2))) & j>2)  
       if(j<257)
            beforind = 1;
            afterind = 1;
            if ~isnan(cell2mat(Data(j,2)))
                temp = cell2mat(Data(j,2));
                while isnan(cell2mat(Data(j+afterind,2))) & j<(j+afterind)-1 
                    afterind=afterind+1;
                end
                while isnan(cell2mat(Data(j-beforind,2))) & j>(j-beforind)+1
                    beforind=beforind+1;
                end
                if (j-beforind)>2 & ~isnan(cell2mat(Data(j-beforind,2)))
                    mat1 = cell2mat(Data(j-beforind,2));
                else
                    beforind = -1;
                    mat1 = cell2mat(Data(j,2));
                end
                if (j+afterind)<=length(Data(:,2)) & ~isnan(cell2mat(Data(j+afterind,2)))
                    mat2 = cell2mat(Data(j+afterind,2));
                else
                    afterind = -1;
                    mat2 = cell2mat(Data(j,2));
                end
                if(afterind==-1)
                    Data(j,11) = num2cell((temp-mat1)/beforind);
                elseif(beforind==-1)
                    Data(j,11) = num2cell((mat2-temp)/afterind);
                else
                    Data(j,11) = num2cell(mean([(mat2-temp)/afterind,(temp-mat1)/beforind]));
                end
                if (j-beforind)>2 & ~isnan(cell2mat(Data(j-beforind,4)))
                    mat1 = cell2mat(Data(j-beforind,4));
                else
                    beforind = -1;
                    mat1 = cell2mat(Data(j,4));
                end
                if (j+afterind)<=length(Data(:,2)) & ~isnan(cell2mat(Data(j+afterind,4)))
                    mat2 = cell2mat(Data(j+afterind,4));
                else
                    afterind = -1;
                    mat2 = cell2mat(Data(j,4));
                end
                temp = cell2mat(Data(j,4));
                if(afterind==-1) 
                    Data(j,12) = num2cell((temp-mat1)/beforind);
                elseif(beforind==-1)
                    Data(j,12) = num2cell((mat2-temp)/afterind);
                else
                    Data(j,12) = num2cell(mean( [(mat2-temp)/afterind,(temp-mat1)/beforind]));
                end
            end
       else
            beforind = 1;
            if ~isnan(cell2mat(Data(j,2)))
                while isnan(cell2mat(Data(j-beforind,2)))
                    beforind=beforind+1;
                end
                if (j-beforind)>2 & ~isnan(cell2mat(Data(j-beforind,2)))
                    mat1 = cell2mat(Data(j-beforind,2));
                else
                     mat1 = cell2mat(Data(j,2));
                end
                mat2 = cell2mat(Data(j,2));
                temp = mat2 - mat1;
                Data(j,11) = num2cell((temp)/(beforind));
                if (j-beforind)>2 & ~isnan(cell2mat(Data(j-beforind,4)))
                    mat1 = cell2mat(Data(j-beforind,4));
                else
                    mat1 = cell2mat(Data(j,4));
                end
                mat2 = cell2mat(Data(j,4));
                temp = mat2 - mat1;
                Data(j,12) = num2cell((temp)/(beforind));
            end
       end
    end
end

% ST3 

temp1 = (ST3(Data(2:end,11),St3_m,St3_P));
temp2 = zeros(1,length(Data(2:end,11))-length(temp1));
temp = [temp2,temp1]';
Data(2:end,13) = num2cell(temp);
temp1 = (ST3(Data(2:end,12),St3_m,St3_P));
temp2 = zeros(1,length(Data(2:end,12))-length(temp1));
temp = [temp2,temp1]';
Data(2:end,14) = num2cell(temp);

% Fitting exponential equations

[R1,R2,X_,Y_] = ExpFitting(Data(3:end,1),Data(3:end,2));
 Data(2,15) = num2cell(R2.rsquare);
 Data(2,16) = num2cell(R2.rmse);
 Data(2,17) = num2cell(R1.a);
 Data(2,18) = num2cell(R1.b);
[R1,R2,X_,Y_] = ExpFitting(Data(3:end,3),Data(3:end,4));
 Data(2,19) = num2cell(R2.rsquare);
 Data(2,20) = num2cell(R2.rmse);
 Data(2,21) = num2cell(R1.a);
 
 Data(2,22) = num2cell(R1.b); 

% Writing single tables

T = table(Data(2:end,1), Data(2:end,2), Data(2:end,3), Data(2:end,4), Data(2:end,5),...
    Data(2:end,6), Data(2:end,7),Data(2:end,8), Data(2:end,9), Data(2:end,10),Data(2:end,11),Data(2:end,12),Data(2:end,13),Data(2:end,14),...
    Data(2:end,15),Data(2:end,16),Data(2:end,17),Data(2:end,18),Data(2:end,19),Data(2:end,20),Data(2:end,21),Data(2:end,22));
T.Properties.VariableNames={'Green_level','Mean_Red','Green','Mean_Blue','Var_Red','Var_Blue','Number_of_pixels','Red_AUC','Blue_AUC','RB_ABC'...
    ,'Red_Slope','Blue_slope','Red_ST3_class','Blue_ST3_class','R2_exp_Red','RMSE_exp_Red ','exp_a_Red','exp_b_Red',...
    'R2_exp_Blue','RMSE_exp_Blue ','exp_a_Blue','exp_b_Blue'};
writetable(T,['./Results/Results of image ',[images(i-1).name],'.csv'],'Delimiter',',');
Data_total_st2(i,1) = cellstr(images(i-1).name);
for ii=2:length(Data(:,2))-1
    Data_total_st2(i,ii)=Data(ii+1,2);
end

for ii=2:length(Data(:,4))-1
    Data_total_st2(i,ii+255)=Data(ii+1,4);
end

Data_total_Exp(i,1) = cellstr(images(i-1).name);
Data_total_Exp(i,2) = Data(2,15);
Data_total_Exp(i,3) = Data(2,16);
Data_total_Exp(i,4) = Data(2,17);
Data_total_Exp(i,5) = Data(2,18);
Data_total_Exp(i,6) = Data(2,19);
Data_total_Exp(i,7) = Data(2,20);
Data_total_Exp(i,8) = Data(2,21);
Data_total_Exp(i,9) = Data(2,22);

GSMRow = zeros(length(Data(3:end,1)),1);
for k=3:length(Data(:,1))
    GSMRow(k,1) = cell2mat(Data(k,1)); 
end

GSMMean_Red = zeros(length(Data(:,2)),1);
for k=3:length(Data(:,2))
    GSMMean_Red(k,1) = cell2mat(Data(k,2)); 
end

GSMGreen = zeros(length(Data(:,3)),1);
for k=3:length(Data(:,3))
    GSMGreen(k,1) = cell2mat(Data(k,3)); 
end

GSMMean_Blue = zeros(length(Data(:,4)),1);
for k=3:length(Data(:,4))
    GSMMean_Blue(k,1) = cell2mat(Data(k,4)); 
end

GSMNO_of_Pixel = [];
for k=3:length(Data(:,7))
    GSMNO_of_Pixel = [GSMNO_of_Pixel;cell2mat(Data(k,7))]; 
end

% Graphs

h = figure();
ax1 = subplot(2,1,1);
hold on
xlim manual
xlim([-1 256])
ylim manual
ylim([-1 256])
plot(GSMRow,GSMMean_Red,'.r');
plot(GSMRow,GSMGreen,'.g');
plot(GSMRow,GSMMean_Blue,'.b');
hold off
title(ax1,'GSM Graph');
ylabel(ax1,'Mean Red, Green, Mean Blue (0-255)');
xlabel(ax1,'Green Level (1-255)');
pbaspect([1 1 1])
ax2 = subplot(2,1,2);
hold on
xlim manual
xlim([-1 256])
ylim auto
plot(GSMRow(3:end,1),GSMNO_of_Pixel,'.k');
title(ax2,'Number of pixels');
ylabel(ax2,'Number of pixels');
xlabel(ax2,'Green Level (1-255)');
pbaspect([1 1 1])
hold off
set(gcf, 'PaperUnits', 'inches');
x_width=7.5; y_width=9.25;
set(gcf, 'PaperPosition', [0 0 x_width y_width]);
saveas(h,['./Graphs/Graphs of image ',[images(i-1).name],'.pdf'],'pdf')
close all
end

T_Exp = table(Data_total_Exp(2:end,1),Data_total_Exp(2:end,2),Data_total_Exp(2:end,3),Data_total_Exp(2:end,4),Data_total_Exp(2:end,5),...
Data_total_Exp(2:end,6),Data_total_Exp(2:end,7),Data_total_Exp(2:end,8),Data_total_Exp(2:end,9));
T_Exp.Properties.VariableNames={'Image','R2_exp_Red','RMSE_exp_Red','Exp_a_Red','Exp_b_Red','R2_exp_Blue','RMSE_exp_Blue','Exp_a_Blue','Exp_b_Blue'};
writetable(T_Exp,'./Results/Total GSM Exp.csv','Delimiter',',');


T_St2 = table(Data_total_st2( 2:end ,1),Data_total_st2( 2:end ,2),Data_total_st2( 2:end ,3),Data_total_st2( 2:end ,4),Data_total_st2( 2:end ,5),Data_total_st2( 2:end ,6),Data_total_st2( 2:end ,7),...
Data_total_st2( 2:end ,8),Data_total_st2( 2:end ,9),Data_total_st2( 2:end ,10),Data_total_st2( 2:end ,11),Data_total_st2( 2:end ,12),Data_total_st2( 2:end ,13),...
Data_total_st2( 2:end ,14),Data_total_st2( 2:end ,15),Data_total_st2( 2:end ,16),Data_total_st2( 2:end ,17),Data_total_st2( 2:end ,18),Data_total_st2( 2:end ,19),...
Data_total_st2( 2:end ,20),Data_total_st2( 2:end ,21),Data_total_st2( 2:end ,22),Data_total_st2( 2:end ,23),Data_total_st2( 2:end ,24),Data_total_st2( 2:end ,25),...
Data_total_st2( 2:end ,26),Data_total_st2( 2:end ,27),Data_total_st2( 2:end ,28),Data_total_st2( 2:end ,29),Data_total_st2( 2:end ,30),Data_total_st2( 2:end ,31),...
Data_total_st2( 2:end ,32),Data_total_st2( 2:end ,33),Data_total_st2( 2:end ,34),Data_total_st2( 2:end ,35),Data_total_st2( 2:end ,36),Data_total_st2( 2:end ,37),...
Data_total_st2( 2:end ,38),Data_total_st2( 2:end ,39),Data_total_st2( 2:end ,40),Data_total_st2( 2:end ,41),Data_total_st2( 2:end ,42),Data_total_st2( 2:end ,43),...
Data_total_st2( 2:end ,44),Data_total_st2( 2:end ,45),Data_total_st2( 2:end ,46),Data_total_st2( 2:end ,47),Data_total_st2( 2:end ,48),Data_total_st2( 2:end ,49),...
Data_total_st2( 2:end ,50),Data_total_st2( 2:end ,51),Data_total_st2( 2:end ,52),Data_total_st2( 2:end ,53),Data_total_st2( 2:end ,54),Data_total_st2( 2:end ,55),...
Data_total_st2( 2:end ,56),Data_total_st2( 2:end ,57),Data_total_st2( 2:end ,58),Data_total_st2( 2:end ,59),Data_total_st2( 2:end ,60),Data_total_st2( 2:end ,61),...
Data_total_st2( 2:end ,62),Data_total_st2( 2:end ,63),Data_total_st2( 2:end ,64),Data_total_st2( 2:end ,65),Data_total_st2( 2:end ,66),Data_total_st2( 2:end ,67),...
Data_total_st2( 2:end ,68),Data_total_st2( 2:end ,69),Data_total_st2( 2:end ,70),Data_total_st2( 2:end ,71),Data_total_st2( 2:end ,72),Data_total_st2( 2:end ,73),...
Data_total_st2( 2:end ,74),Data_total_st2( 2:end ,75),Data_total_st2( 2:end ,76),Data_total_st2( 2:end ,77),Data_total_st2( 2:end ,78),Data_total_st2( 2:end ,79),...
Data_total_st2( 2:end ,80),Data_total_st2( 2:end ,81),Data_total_st2( 2:end ,82),Data_total_st2( 2:end ,83),Data_total_st2( 2:end ,84),Data_total_st2( 2:end ,85),...
Data_total_st2( 2:end ,86),Data_total_st2( 2:end ,87),Data_total_st2( 2:end ,88),Data_total_st2( 2:end ,89),Data_total_st2( 2:end ,90),Data_total_st2( 2:end ,91),...
Data_total_st2( 2:end ,92),Data_total_st2( 2:end ,93),Data_total_st2( 2:end ,94),Data_total_st2( 2:end ,95),Data_total_st2( 2:end ,96),Data_total_st2( 2:end ,97),...
Data_total_st2( 2:end ,98),Data_total_st2( 2:end ,99),Data_total_st2( 2:end ,100),Data_total_st2( 2:end ,101),Data_total_st2( 2:end ,102),Data_total_st2( 2:end ,103),...
Data_total_st2( 2:end ,104),Data_total_st2( 2:end ,105),Data_total_st2( 2:end ,106),Data_total_st2( 2:end ,107),Data_total_st2( 2:end ,108),Data_total_st2( 2:end ,109),...
Data_total_st2( 2:end ,110),Data_total_st2( 2:end ,111),Data_total_st2( 2:end ,112),Data_total_st2( 2:end ,113),Data_total_st2( 2:end ,114),Data_total_st2( 2:end ,115),...
Data_total_st2( 2:end ,116),Data_total_st2( 2:end ,117),Data_total_st2( 2:end ,118),Data_total_st2( 2:end ,119),Data_total_st2( 2:end ,120),Data_total_st2( 2:end ,121),...
Data_total_st2( 2:end ,122),Data_total_st2( 2:end ,123),Data_total_st2( 2:end ,124),Data_total_st2( 2:end ,125),Data_total_st2( 2:end ,126),Data_total_st2( 2:end ,127),...
Data_total_st2( 2:end ,128),Data_total_st2( 2:end ,129),Data_total_st2( 2:end ,130),Data_total_st2( 2:end ,131),Data_total_st2( 2:end ,132),Data_total_st2( 2:end ,133),...
Data_total_st2( 2:end ,134),Data_total_st2( 2:end ,135),Data_total_st2( 2:end ,136),Data_total_st2( 2:end ,137),Data_total_st2( 2:end ,138),Data_total_st2( 2:end ,139),...
Data_total_st2( 2:end ,140),Data_total_st2( 2:end ,141),Data_total_st2( 2:end ,142),Data_total_st2( 2:end ,143),Data_total_st2( 2:end ,144),Data_total_st2( 2:end ,145),...
Data_total_st2( 2:end ,146),Data_total_st2( 2:end ,147),Data_total_st2( 2:end ,148),Data_total_st2( 2:end ,149),Data_total_st2( 2:end ,150),Data_total_st2( 2:end ,151),...
Data_total_st2( 2:end ,152),Data_total_st2( 2:end ,153),Data_total_st2( 2:end ,154),Data_total_st2( 2:end ,155),Data_total_st2( 2:end ,156),Data_total_st2( 2:end ,157),...
Data_total_st2( 2:end ,158),Data_total_st2( 2:end ,159),Data_total_st2( 2:end ,160),Data_total_st2( 2:end ,161),Data_total_st2( 2:end ,162),Data_total_st2( 2:end ,163),...
Data_total_st2( 2:end ,164),Data_total_st2( 2:end ,165),Data_total_st2( 2:end ,166),Data_total_st2( 2:end ,167),Data_total_st2( 2:end ,168),Data_total_st2( 2:end ,169),...
Data_total_st2( 2:end ,170),Data_total_st2( 2:end ,171),Data_total_st2( 2:end ,172),Data_total_st2( 2:end ,173),Data_total_st2( 2:end ,174),Data_total_st2( 2:end ,175),...
Data_total_st2( 2:end ,176),Data_total_st2( 2:end ,177),Data_total_st2( 2:end ,178),Data_total_st2( 2:end ,179),Data_total_st2( 2:end ,180),Data_total_st2( 2:end ,181),...
Data_total_st2( 2:end ,182),Data_total_st2( 2:end ,183),Data_total_st2( 2:end ,184),Data_total_st2( 2:end ,185),Data_total_st2( 2:end ,186),Data_total_st2( 2:end ,187),...
Data_total_st2( 2:end ,188),Data_total_st2( 2:end ,189),Data_total_st2( 2:end ,190),Data_total_st2( 2:end ,191),Data_total_st2( 2:end ,192),Data_total_st2( 2:end ,193),...
Data_total_st2( 2:end ,194),Data_total_st2( 2:end ,195),Data_total_st2( 2:end ,196),Data_total_st2( 2:end ,197),Data_total_st2( 2:end ,198),Data_total_st2( 2:end ,199),...
Data_total_st2( 2:end ,200),Data_total_st2( 2:end ,201),Data_total_st2( 2:end ,202),Data_total_st2( 2:end ,203),Data_total_st2( 2:end ,204),Data_total_st2( 2:end ,205),...
Data_total_st2( 2:end ,206),Data_total_st2( 2:end ,207),Data_total_st2( 2:end ,208),Data_total_st2( 2:end ,209),Data_total_st2( 2:end ,210),Data_total_st2( 2:end ,211),...
Data_total_st2( 2:end ,212),Data_total_st2( 2:end ,213),Data_total_st2( 2:end ,214),Data_total_st2( 2:end ,215),Data_total_st2( 2:end ,216),Data_total_st2( 2:end ,217),...
Data_total_st2( 2:end ,218),Data_total_st2( 2:end ,219),Data_total_st2( 2:end ,220),Data_total_st2( 2:end ,221),Data_total_st2( 2:end ,222),Data_total_st2( 2:end ,223),...
Data_total_st2( 2:end ,224),Data_total_st2( 2:end ,225),Data_total_st2( 2:end ,226),Data_total_st2( 2:end ,227),Data_total_st2( 2:end ,228),Data_total_st2( 2:end ,229),...
Data_total_st2( 2:end ,230),Data_total_st2( 2:end ,231),Data_total_st2( 2:end ,232),Data_total_st2( 2:end ,233),Data_total_st2( 2:end ,234),Data_total_st2( 2:end ,235),...
Data_total_st2( 2:end ,236),Data_total_st2( 2:end ,237),Data_total_st2( 2:end ,238),Data_total_st2( 2:end ,239),Data_total_st2( 2:end ,240),Data_total_st2( 2:end ,241),...
Data_total_st2( 2:end ,242),Data_total_st2( 2:end ,243),Data_total_st2( 2:end ,244),Data_total_st2( 2:end ,245),Data_total_st2( 2:end ,246),Data_total_st2( 2:end ,247),...
Data_total_st2( 2:end ,248),Data_total_st2( 2:end ,249),Data_total_st2( 2:end ,250),Data_total_st2( 2:end ,251),Data_total_st2( 2:end ,252),Data_total_st2( 2:end ,253),...
Data_total_st2( 2:end ,254),Data_total_st2( 2:end ,255),Data_total_st2( 2:end ,256),Data_total_st2( 2:end ,257),Data_total_st2( 2:end ,258),Data_total_st2( 2:end ,259),...
Data_total_st2( 2:end ,260),Data_total_st2( 2:end ,261),Data_total_st2( 2:end ,262),Data_total_st2( 2:end ,263),Data_total_st2( 2:end ,264),Data_total_st2( 2:end ,265),...
Data_total_st2( 2:end ,266),Data_total_st2( 2:end ,267),Data_total_st2( 2:end ,268),Data_total_st2( 2:end ,269),Data_total_st2( 2:end ,270),Data_total_st2( 2:end ,271),...
Data_total_st2( 2:end ,272),Data_total_st2( 2:end ,273),Data_total_st2( 2:end ,274),Data_total_st2( 2:end ,275),Data_total_st2( 2:end ,276),Data_total_st2( 2:end ,277),...
Data_total_st2( 2:end ,278),Data_total_st2( 2:end ,279),Data_total_st2( 2:end ,280),Data_total_st2( 2:end ,281),Data_total_st2( 2:end ,282),Data_total_st2( 2:end ,283),...
Data_total_st2( 2:end ,284),Data_total_st2( 2:end ,285),Data_total_st2( 2:end ,286),Data_total_st2( 2:end ,287),Data_total_st2( 2:end ,288),Data_total_st2( 2:end ,289),...
Data_total_st2( 2:end ,290),Data_total_st2( 2:end ,291),Data_total_st2( 2:end ,292),Data_total_st2( 2:end ,293),Data_total_st2( 2:end ,294),Data_total_st2( 2:end ,295),...
Data_total_st2( 2:end ,296),Data_total_st2( 2:end ,297),Data_total_st2( 2:end ,298),Data_total_st2( 2:end ,299),Data_total_st2( 2:end ,300),Data_total_st2( 2:end ,301),...
Data_total_st2( 2:end ,302),Data_total_st2( 2:end ,303),Data_total_st2( 2:end ,304),Data_total_st2( 2:end ,305),Data_total_st2( 2:end ,306),Data_total_st2( 2:end ,307),...
Data_total_st2( 2:end ,308),Data_total_st2( 2:end ,309),Data_total_st2( 2:end ,310),Data_total_st2( 2:end ,311),Data_total_st2( 2:end ,312),Data_total_st2( 2:end ,313),...
Data_total_st2( 2:end ,314),Data_total_st2( 2:end ,315),Data_total_st2( 2:end ,316),Data_total_st2( 2:end ,317),Data_total_st2( 2:end ,318),Data_total_st2( 2:end ,319),...
Data_total_st2( 2:end ,320),Data_total_st2( 2:end ,321),Data_total_st2( 2:end ,322),Data_total_st2( 2:end ,323),Data_total_st2( 2:end ,324),Data_total_st2( 2:end ,325),...
Data_total_st2( 2:end ,326),Data_total_st2( 2:end ,327),Data_total_st2( 2:end ,328),Data_total_st2( 2:end ,329),Data_total_st2( 2:end ,330),Data_total_st2( 2:end ,331),...
Data_total_st2( 2:end ,332),Data_total_st2( 2:end ,333),Data_total_st2( 2:end ,334),Data_total_st2( 2:end ,335),Data_total_st2( 2:end ,336),Data_total_st2( 2:end ,337),...
Data_total_st2( 2:end ,338),Data_total_st2( 2:end ,339),Data_total_st2( 2:end ,340),Data_total_st2( 2:end ,341),Data_total_st2( 2:end ,342),Data_total_st2( 2:end ,343),...
Data_total_st2( 2:end ,344),Data_total_st2( 2:end ,345),Data_total_st2( 2:end ,346),Data_total_st2( 2:end ,347),Data_total_st2( 2:end ,348),Data_total_st2( 2:end ,349),...
Data_total_st2( 2:end ,350),Data_total_st2( 2:end ,351),Data_total_st2( 2:end ,352),Data_total_st2( 2:end ,353),Data_total_st2( 2:end ,354),Data_total_st2( 2:end ,355),...
Data_total_st2( 2:end ,356),Data_total_st2( 2:end ,357),Data_total_st2( 2:end ,358),Data_total_st2( 2:end ,359),Data_total_st2( 2:end ,360),Data_total_st2( 2:end ,361),...
Data_total_st2( 2:end ,362),Data_total_st2( 2:end ,363),Data_total_st2( 2:end ,364),Data_total_st2( 2:end ,365),Data_total_st2( 2:end ,366),Data_total_st2( 2:end ,367),...
Data_total_st2( 2:end ,368),Data_total_st2( 2:end ,369),Data_total_st2( 2:end ,370),Data_total_st2( 2:end ,371),Data_total_st2( 2:end ,372),Data_total_st2( 2:end ,373),...
Data_total_st2( 2:end ,374),Data_total_st2( 2:end ,375),Data_total_st2( 2:end ,376),Data_total_st2( 2:end ,377),Data_total_st2( 2:end ,378),Data_total_st2( 2:end ,379),...
Data_total_st2( 2:end ,380),Data_total_st2( 2:end ,381),Data_total_st2( 2:end ,382),Data_total_st2( 2:end ,383),Data_total_st2( 2:end ,384),Data_total_st2( 2:end ,385),...
Data_total_st2( 2:end ,386),Data_total_st2( 2:end ,387),Data_total_st2( 2:end ,388),Data_total_st2( 2:end ,389),Data_total_st2( 2:end ,390),Data_total_st2( 2:end ,391),...
Data_total_st2( 2:end ,392),Data_total_st2( 2:end ,393),Data_total_st2( 2:end ,394),Data_total_st2( 2:end ,395),Data_total_st2( 2:end ,396),Data_total_st2( 2:end ,397),...
Data_total_st2( 2:end ,398),Data_total_st2( 2:end ,399),Data_total_st2( 2:end ,400),Data_total_st2( 2:end ,401),Data_total_st2( 2:end ,402),Data_total_st2( 2:end ,403),...
Data_total_st2( 2:end ,404),Data_total_st2( 2:end ,405),Data_total_st2( 2:end ,406),Data_total_st2( 2:end ,407),Data_total_st2( 2:end ,408),Data_total_st2( 2:end ,409),...
Data_total_st2( 2:end ,410),Data_total_st2( 2:end ,411),Data_total_st2( 2:end ,412),Data_total_st2( 2:end ,413),Data_total_st2( 2:end ,414),Data_total_st2( 2:end ,415),...
Data_total_st2( 2:end ,416),Data_total_st2( 2:end ,417),Data_total_st2( 2:end ,418),Data_total_st2( 2:end ,419),Data_total_st2( 2:end ,420),Data_total_st2( 2:end ,421),...
Data_total_st2( 2:end ,422),Data_total_st2( 2:end ,423),Data_total_st2( 2:end ,424),Data_total_st2( 2:end ,425),Data_total_st2( 2:end ,426),Data_total_st2( 2:end ,427),...
Data_total_st2( 2:end ,428),Data_total_st2( 2:end ,429),Data_total_st2( 2:end ,430),Data_total_st2( 2:end ,431),Data_total_st2( 2:end ,432),Data_total_st2( 2:end ,433),...
Data_total_st2( 2:end ,434),Data_total_st2( 2:end ,435),Data_total_st2( 2:end ,436),Data_total_st2( 2:end ,437),Data_total_st2( 2:end ,438),Data_total_st2( 2:end ,439),...
Data_total_st2( 2:end ,440),Data_total_st2( 2:end ,441),Data_total_st2( 2:end ,442),Data_total_st2( 2:end ,443),Data_total_st2( 2:end ,444),Data_total_st2( 2:end ,445),...
Data_total_st2( 2:end ,446),Data_total_st2( 2:end ,447),Data_total_st2( 2:end ,448),Data_total_st2( 2:end ,449),Data_total_st2( 2:end ,450),Data_total_st2( 2:end ,451),...
Data_total_st2( 2:end ,452),Data_total_st2( 2:end ,453),Data_total_st2( 2:end ,454),Data_total_st2( 2:end ,455),Data_total_st2( 2:end ,456),Data_total_st2( 2:end ,457),...
Data_total_st2( 2:end ,458),Data_total_st2( 2:end ,459),Data_total_st2( 2:end ,460),Data_total_st2( 2:end ,461),Data_total_st2( 2:end ,462),Data_total_st2( 2:end ,463),...
Data_total_st2( 2:end ,464),Data_total_st2( 2:end ,465),Data_total_st2( 2:end ,466),Data_total_st2( 2:end ,467),Data_total_st2( 2:end ,468),Data_total_st2( 2:end ,469),...
Data_total_st2( 2:end ,470),Data_total_st2( 2:end ,471),Data_total_st2( 2:end ,472),Data_total_st2( 2:end ,473),Data_total_st2( 2:end ,474),Data_total_st2( 2:end ,475),...
Data_total_st2( 2:end ,476),Data_total_st2( 2:end ,477),Data_total_st2( 2:end ,478),Data_total_st2( 2:end ,479),Data_total_st2( 2:end ,480),Data_total_st2( 2:end ,481),...
Data_total_st2( 2:end ,482),Data_total_st2( 2:end ,483),Data_total_st2( 2:end ,484),Data_total_st2( 2:end ,485),Data_total_st2( 2:end ,486),Data_total_st2( 2:end ,487),...
Data_total_st2( 2:end ,488),Data_total_st2( 2:end ,489),Data_total_st2( 2:end ,490),Data_total_st2( 2:end ,491),Data_total_st2( 2:end ,492),Data_total_st2( 2:end ,493),...
Data_total_st2( 2:end ,494),Data_total_st2( 2:end ,495),Data_total_st2( 2:end ,496),Data_total_st2( 2:end ,497),Data_total_st2( 2:end ,498),Data_total_st2( 2:end ,499),...
Data_total_st2( 2:end ,500),Data_total_st2( 2:end ,501),Data_total_st2( 2:end ,502),Data_total_st2( 2:end ,503),Data_total_st2( 2:end ,504),Data_total_st2( 2:end ,505),...
Data_total_st2( 2:end ,506),Data_total_st2( 2:end ,507),Data_total_st2( 2:end ,508),Data_total_st2( 2:end ,509),Data_total_st2( 2:end ,510),Data_total_st2( 2:end ,511));
T_St2.Properties.VariableNames={'Image','R_G1','R_G2','R_G3','R_G4','R_G5','R_G6','R_G7','R_G8','R_G9','R_G10','R_G11','R_G12','R_G13',...
'R_G14','R_G15','R_G16','R_G17','R_G18','R_G19','R_G20','R_G21','R_G22','R_G23','R_G24','R_G25',...
'R_G26','R_G27','R_G28','R_G29','R_G30','R_G31','R_G32','R_G33','R_G34','R_G35','R_G36','R_G37',...
'R_G38','R_G39','R_G40','R_G41','R_G42','R_G43','R_G44','R_G45','R_G46','R_G47','R_G48','R_G49',...
'R_G50','R_G51','R_G52','R_G53','R_G54','R_G55','R_G56','R_G57','R_G58','R_G59','R_G60','R_G61',...
'R_G62','R_G63','R_G64','R_G65','R_G66','R_G67','R_G68','R_G69','R_G70','R_G71','R_G72','R_G73',...
'R_G74','R_G75','R_G76','R_G77','R_G78','R_G79','R_G80','R_G81','R_G82','R_G83','R_G84','R_G85',...
'R_G86','R_G87','R_G88','R_G89','R_G90','R_G91','R_G92','R_G93','R_G94','R_G95','R_G96','R_G97',...
'R_G98','R_G99','R_G100','R_G101','R_G102','R_G103','R_G104','R_G105','R_G106','R_G107','R_G108',...
'R_G109','R_G110','R_G111','R_G112','R_G113','R_G114','R_G115','R_G116','R_G117','R_G118','R_G119',...
'R_G120','R_G121','R_G122','R_G123','R_G124','R_G125','R_G126','R_G127','R_G128','R_G129','R_G130',...
'R_G131','R_G132','R_G133','R_G134','R_G135','R_G136','R_G137','R_G138','R_G139','R_G140','R_G141',...
'R_G142','R_G143','R_G144','R_G145','R_G146','R_G147','R_G148','R_G149','R_G150','R_G151','R_G152',...
'R_G153','R_G154','R_G155','R_G156','R_G157','R_G158','R_G159','R_G160','R_G161','R_G162','R_G163',...
'R_G164','R_G165','R_G166','R_G167','R_G168','R_G169','R_G170','R_G171','R_G172','R_G173','R_G174',...
'R_G175','R_G176','R_G177','R_G178','R_G179','R_G180','R_G181','R_G182','R_G183','R_G184','R_G185',...
'R_G186','R_G187','R_G188','R_G189','R_G190','R_G191','R_G192','R_G193','R_G194','R_G195','R_G196',...
'R_G197','R_G198','R_G199','R_G200','R_G201','R_G202','R_G203','R_G204','R_G205','R_G206','R_G207',...
'R_G208','R_G209','R_G210','R_G211','R_G212','R_G213','R_G214','R_G215','R_G216','R_G217','R_G218',...
'R_G219','R_G220','R_G221','R_G222','R_G223','R_G224','R_G225','R_G226','R_G227','R_G228','R_G229',...
'R_G230','R_G231','R_G232','R_G233','R_G234','R_G235','R_G236','R_G237','R_G238','R_G239','R_G240',...
'R_G241','R_G242','R_G243','R_G244','R_G245','R_G246','R_G247','R_G248','R_G249','R_G250','R_G251',...
'R_G252','R_G253','R_G254','R_G255','B_G1','B_G2','B_G3','B_G4','B_G5','B_G6','B_G7','B_G8','B_G9',...
'B_G10','B_G11','B_G12','B_G13','B_G14','B_G15','B_G16','B_G17','B_G18','B_G19','B_G20','B_G21',...
'B_G22','B_G23','B_G24','B_G25','B_G26','B_G27','B_G28','B_G29','B_G30','B_G31','B_G32','B_G33',...
'B_G34','B_G35','B_G36','B_G37','B_G38','B_G39','B_G40','B_G41','B_G42','B_G43','B_G44','B_G45',...
'B_G46','B_G47','B_G48','B_G49','B_G50','B_G51','B_G52','B_G53','B_G54','B_G55','B_G56','B_G57',...
'B_G58','B_G59','B_G60','B_G61','B_G62','B_G63','B_G64','B_G65','B_G66','B_G67','B_G68','B_G69',...
'B_G70','B_G71','B_G72','B_G73','B_G74','B_G75','B_G76','B_G77','B_G78','B_G79','B_G80','B_G81',...
'B_G82','B_G83','B_G84','B_G85','B_G86','B_G87','B_G88','B_G89','B_G90','B_G91','B_G92','B_G93',...
'B_G94','B_G95','B_G96','B_G97','B_G98','B_G99','B_G100','B_G101','B_G102','B_G103','B_G104',...
'B_G105','B_G106','B_G107','B_G108','B_G109','B_G110','B_G111','B_G112','B_G113','B_G114','B_G115',...
'B_G116','B_G117','B_G118','B_G119','B_G120','B_G121','B_G122','B_G123','B_G124','B_G125','B_G126',...
'B_G127','B_G128','B_G129','B_G130','B_G131','B_G132','B_G133','B_G134','B_G135','B_G136','B_G137',...
'B_G138','B_G139','B_G140','B_G141','B_G142','B_G143','B_G144','B_G145','B_G146','B_G147','B_G148',...
'B_G149','B_G150','B_G151','B_G152','B_G153','B_G154','B_G155','B_G156','B_G157','B_G158','B_G159',...
'B_G160','B_G161','B_G162','B_G163','B_G164','B_G165','B_G166','B_G167','B_G168','B_G169','B_G170',...
'B_G171','B_G172','B_G173','B_G174','B_G175','B_G176','B_G177','B_G178','B_G179','B_G180','B_G181',...
'B_G182','B_G183','B_G184','B_G185','B_G186','B_G187','B_G188','B_G189','B_G190','B_G191','B_G192',...
'B_G193','B_G194','B_G195','B_G196','B_G197','B_G198','B_G199','B_G200','B_G201','B_G202','B_G203',...
'B_G204','B_G205','B_G206','B_G207','B_G208','B_G209','B_G210','B_G211','B_G212','B_G213','B_G214',...
'B_G215','B_G216','B_G217','B_G218','B_G219','B_G220','B_G221','B_G222','B_G223','B_G224','B_G225',...
'B_G226','B_G227','B_G228','B_G229','B_G230','B_G231','B_G232','B_G233','B_G234','B_G235','B_G236',...
'B_G237','B_G238','B_G239','B_G240','B_G241','B_G242','B_G243','B_G244','B_G245','B_G246','B_G247',...
'B_G248','B_G249','B_G250','B_G251','B_G252','B_G253','B_G254','B_G255'};
writetable(T_St2,'./Results/Total GSM ST2.csv','Delimiter',',');


disp('===== Processing completed! =====')

% End