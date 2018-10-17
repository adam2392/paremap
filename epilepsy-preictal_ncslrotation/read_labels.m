%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% m-file: read_labels.m
%
% Description: It stores in a struct array the following info for each
%              subject in the repository:
%               - id of the subject;
%               - labels of the channels used for SVD;
%               - labels of the focal channels;
%               - position of the channels used for SVD in the files of raw
%                 recordings.
%
%              The struct array is stored in a .mat file. Note that, for
%              those subjects for which the number of raw channels varies
%              across the days, only the map of the channels during the day
%              of seizures is stored.
%
%
%
% Author: S. Santaniello
% Modified by: B. Chennuri
%   
% Ver.: 4.0 - Date: 02/03/2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
clc

if exist('infolabels.mat','file')
    delete('infolabels.mat');
end

% fill in the struct array
labels(1).subject  = 'PY04N007';
labels(1).values   = {'LFG1','LFG2','LTG41','LTG42','LFG3','LFG4','LFG5','LFG6','LFG7','LFG8','LFG9','LFG10','LFG11','LFG12','LFG13','LFG14','LFG15','LFG16','LFG17','LFG18','LFG19','LFG20','LTG1','LTG2','LTG3','LTG4','LTG5','LTG6','LTG7','LTG8','LTG9','LTG10','LTG11','LTG12','LTG13','LTG14','LTG15','LTG16','LTG17','LTG18','LTG19','LTG20','LTG21','LTG22','LTG23','LTG24','LTG25','LTG26','LTG27','LTG28','LTG29','LTG30','LTG31','LTG32','LTG33','LTG34','LTG35','LTG36','LTG37','LTG38','LTG39','LTG40','LTG45','LTG46','LTG47','LTG48','LIT1','LIT2','LIT3','LIT4','LIT5','LIT6','LIT7','LIT8','LIT9','LIT10','LAT1','LAT2','LAT3','LMT1','LMT2','LMT3','LMT4','LPT1','LPT2','LPT3'};
labels(1).focus    = {'LIT1','LIT2','LIT6','LIT7','LAT1','LAT2','LAT3'};
labels(1).channel  = [1:79 81:87];


labels(2).subject  = 'PY04N008';
labels(2).values   = {'IPO1','IPO2','IPO3','IPO4','IPO5','IPO6','IPO7','IPO8','IPO9','IPO10','IPO11','IPO12','IPO13','IPO14','IPO15','IPO16','SPO1','SPO2','SPO3','SPO4','SPO5','SPO6','SPO7','SPO8','ROG1','ROG2','ROG3','ROG4','ROG5','ROG6','ROG7','ROG8','ROG9','ROG10','ROG11','ROG12','ROG13','ROG14','ROG15','ROG16'};
labels(2).focus    = {'IPO1','IPO2','IPO3','IPO9','IPO10','IPO11'};
labels(2).channel  = 1:40;


labels(3).subject  = 'PY04N012';
labels(3).values   = {'LPS1','LPS2','LLF1','LLF2','LPS3','LPS4','LPS5','LPS6','LPS7','LPS8','LPS9','LPS10','LPS11','LPS12','LPS13','LPS14','LPS15','LPS16','LPS17','LPS18','LPS19','LPS20','LPS21','LPS22','LPS23','LPS24','LPS25','LPS26','LPS27','LPS28','LPS29','LPS30','LPS31','LPS32','LPT1','LPT2','LPT3','LPT4','LPT5','LPT6','LPT7','LPT8','LMT1','LMT2','LMT3','LMT4','LMT5','LMT6','LMT7','LMT8','LMT9','LMT10','LAT1','LAT2','LAT3','LAT4','LOF1','LOF2','LOF3','LOF4','LOF5','LOF6','LOF7','LOF8','LOF9','LOF10','LLF5','LLF6','LLF7','LLF8','LLF9','LLF10','LPF1','LPF2','LPF3','LPF4','LPF5','LPF6','LPF7','LPF8','LPF9','LPF10'};
labels(3).focus    = {'LMT1','LMT6'};
labels(3).channel  = 1:82;

labels(4).subject  = 'PY04N013';
labels(4).values   = {'RTG1','RTG2','RTG41','RTG42','RTG3','RTG4','RTG5','RTG6','RTG7','RTG8','RTG10','RTG11','RTG12','RTG13','RTG14','RTG15','RTG16','RTG17','RTG18','RTG19','RTG20','RTG21','RTG22','RTG23','RTG24','RTG25','RTG26','RTG27','RTG28','RTG29','RTG30','RTG31','RTG32','RTG33','RTG34','RTG35','RTG36','RTG37','RTG38','RTG39','RTG40','RTG44','RTG45','RTG46','RTG47','RTG48','RMO1','RMO2','RMO3','RMO4','RMO5','RMO6','RMO7','RMO8','RMO9','RMO10','RPO1','RPO2','RPO3','RPO4','RPO6','RPO8','TPO2','TPO3','TPO4','ROG1','ROG2','ROG3','ROG4','ROG5','ROG6','ROG7','ROG8','ROG9','ROG10','ROG11','ROG12','ROG13','ROG14','ROG15','ROG16','ROG17','ROG18'};
labels(4).focus    = {'RTG3','RTG4','RTG5','RTG6','RTG10','RTG11','RTG12','RTG13','RTG14','RTG15','RTG18','RTG19','RTG20','RTG21','RTG22','RMO1','RMO2','RMO3','RMO4','RMO5','RMO6','RMO7','RMO8','RMO9','RMO10','RPO1','RPO2','RPO3','RPO4','RPO6','RPO8','TPO2','TPO3','TPO4','ROG1','ROG2','ROG3','ROG4','ROG5','ROG6','ROG7','ROG8','ROG9','ROG10'};
labels(4).channel  = [1:10 12:61 63 65 69:89];

labels(5).subject  = 'PY04N015';
labels(5).values   = {'LTG01','LTG02','LTG03','LTG04','LTG05','LTG06','LTG09','LTG10','LTG11','LTG12','LTG13','LTG14','LTG15','LTG16','LTG17','LTG18','LTG19','LTG20','LTG21','LTG22','LTG23','LTG24','LTG25','LTG26','LTG27','LTG28','LTG29','LTG30','LTG31','LTG32','LTG33','LTG34','LTG35','LTG36','LTG37','LTG38','LTG39','LTG40','LTG41','LTG42','LTG43','LTG44','LTG45','LTG46','LTG47','LTG48','LAT01','LAT02','LAT03','LAT04','LAT05','LAT06','LAT07','LAT08','LAT09','LAT10','LAT11','LAT12','LAT13','LAT14','LAT15','LAT16','LPT01','LPT02','LPT03','LPT04','LPT05','LPT06','LPT07','LPT08','LPT09','LPT10','LOF01','LOF02','LOF03','LOF04','LOF05','LOF06','LOF07','LOF08','LOF09','RPT01','RPT02','RPT03','RPT04','RPT05','RAT01','RAT02','RAT03','RAT04','RAT05','RAT06'};
labels(5).focus    = {'LAT01','LAT02','LAT03','LAT04','LAT05','LAT06','LAT07','LAT08','LAT09','LAT10','LAT11','LAT12','LAT13','LAT14','LAT15','LAT16','LPT01','LPT02','LPT03','LPT04','LPT05','LPT06','LPT07','LPT08','LPT09','LPT10'};
labels(5).channel  = [1:2 5:83 85:95];

labels(6).subject  = 'PY05N004';
labels(6).values   = {'LPG01','LPG02','LPG03','LPG04','LPG05','LPG06','LPG07','LPG08','LPG09','LPG10','LPG11','LPG12','LPG14','LPG17','LPG18','LPG19','LPG20','LPG21','LPG22','LPG23','LPG24','LPG25','LPG26','LPG27','LPG28','LPG29','LPG30','LPG31','LPG32','LPG33','LPG34','LPG35','LPG36','LPG37','LPG39','LPG40','LPG41','LPG42','LPG43','LPG44','LPG45','LPG46','LPG47','LPG48','LPG49','LPG50','LPG51','LPG52','LPG53','LPG54','LPG55','LPG56','LPG57','LPG58','LPG59','LPG60','LPG61','LPG62','LPG63','LFG01','LFG02','LFG03','LFG04','LFG05','LFG06','LFG07','LFG08','LFG09','LFG10','LFG11','LFG12','LFG13','LFG14','LFG15','LFG16','LFG17','LFG18','LFG20','LOG01','LOG02','LOG03','LOG04','LOG05','LOG06','LOG07','LOG08','LOG09','LOG10','LOG11','LOG12','LOG13','LOG14','LOG15','LOG16'};
labels(6).focus    = {'LOG02','LOG03','LOG04','LOG05','LOG06','LOG07','LOG10','LOG11','LOG12','LOG13','LOG14','LOG15'};
labels(6).channel  = [1:12 14 17:37 39:81 83:99];

labels(7).subject  = 'PY05N005';
labels(7).values   = {'LAT01','LAT02','LAT03','LAT04','LAT05','LAT06','LAT07','LAT08','LAT09','LAT10','LPT01','LPT02','LPT03','LPT04','LPT05','LPT06','LPT07','LPT08','LPT09','LPT10','LOG01','LOG03','LOG04','LOG05','LOG06','LOG07','LOG08','LOG09','LOG10','LOG11','LOG12','LOG13','LOG14','LOG15','LOG16','LOG17','LOG18','LPO01','LPO02','LPO03','LPO04','LPO05','LPO06','LPO07','LPO08','LPO09','LPO10','LFG31','LFG32','LFG33','LFG34','LFG35','LFG36','LFG37','LFG38','LFG39','LFG40','LFG41','LFG42','LFG43','LFG44','LFG45','LFG46','LFG47','LFG48','LFG49','LFG50','LFG51','LFG52','LFG53','LFG54','LFG55','LFG56','LFG57','LFG58','LFG59','LFG60','LFG61','LFG62','LFG63','LFG01','LFG02','LFG03','LFG04','LFG05','LFG06','LFG07','LFG08','LFG09','LFG10','LFG11','LFG12','LFG13','LFG14','LFG15','LFG16','LFG17','LFG18','LFG19','LFG20','LFG21','LFG22','LFG23','LFG24','LFG25','LFG26','LFG27','LFG28','LFG29','LFG30'};
labels(7).focus    = {'LAT01','LAT02','LAT03','LAT04','LAT05','LAT06','LAT07','LAT08','LAT09','LAT10','LPT06','LPT07','LPT08','LPT09','LPT10','LFG31','LFG32','LFG02','LFG03','LFG04','LFG05','LFG06','LFG07','LFG08','LFG10','LFG11','LFG12','LFG13','LFG14','LFG15','LFG16','LFG18','LFG19','LFG20','LFG21','LFG22','LFG23','LFG24','LFG26','LFG27','LFG28','LFG29','LFG30'};
labels(7).channel  = [1:21 23:38 41:83 89:118];

labels(8).subject  = 'PY11N003';
labels(8).values   = {'LFG3','LFG4','LFG5','LFG6','LFG7','LFG8','LFG9','LFG10','LFG11','LFG12','LFG13','LFG14','LFG15','LFG16','LFG17','LFG18','LFG19','LFG20','LFG21','LFG22','LFG23','LFG24','LFG25','LFG26','LFG27','LFG28','LFG29','LFG30','LFG31','LFG32','LFG33','LFG34','LFG35','LFG36','LFG37','LFG38','LFG39','LFG40','LFG41','LFG42','LFG43','LFG44','LFG45','LFG46','LFG47','LFG48','LFG49','LFG50','LFG51','LFG52','LFG53','LFG54','LFG55','LFG56','LFG57','LFG58','LFG59','LFG60','LFG61','LFG62','LFG63','LFG64','IHG2','IHG3','IHG4','IHG5','IHG6','IHG7','IHG8','IHG9','IHG10','IHG11','IHG12','IHG13','IHG14','IHG15','IHG16','IHG17','IHG18','IHG19','IHG20','IHG21','IHG22','IHG23','IHG24','IHG25','IHG26','IHG27','IHG28','IHG29','IHG30','IHG31','IHG32','ACD1','ACD2','ACD3','ACD4','ACD5','ACD6','ACD7','ACD8','MCD1','MCD2','MCD3','MCD4','MCD5','MCD6','MCD7','MCD8','PCD1','PCD2','PCD3','PCD4','PCD5','PCD6','PCD7','PCD8'};
labels(8).focus    = {'LFG5','LFG6','LFG7','LFG8','LFG13','LFG14','LFG15','LFG16','LFG21','LFG22','LFG23','LFG24','LFG29','LFG30','LFG31','LFG32','LFG37','LFG38','LFG39','LFG40','LFG45','LFG46','LFG47','LFG48','LFG56','LFG64','IHG17','IHG18','IHG19','IHG20','IHG21','IHG22','IHG25','IHG26','IHG27','IHG28','IHG29','IHG30'};
labels(8).channel  = 1:117;

labels(9).subject  = 'PY11N004';
labels(9).values   = {'LFG3','LFG4','LFG5','LFG6','LFG7','LFG9','LFG10','LFG11','LFG12','LFG13','LFG14','LFG15','LFG17','LFG18','LFG19','LFG20','LFG21','LFG22','LFG23','LFG24','LFG25','LFG26','LFG27','LFG28','LFG29','LFG30','LFG31','LFG32','LFG34','LFG35','LFG36','LFG37','LFG38','LFG39','LFG40','LFG41','LFG42','LFG43','LFG44','LFG45','LFG46','LFG47','LFG48','LFG49','LFG50','LFG51','LFG52','LFG53','LFG54','LFG55','LFG56','LFG57','LFG58','LFG59','LFG60','LFG61','LFG62','LFG63','LFG64','IHG1','IHG2','IHG3','IHG4','IHG5','IHG6','IHG7','IHG8','IHG9','IHG10','IHG11','IHG12','IHG13','IHG14','IHG15','IHG16','IHG17','IHG18','IHG19','IHG20','IHG21','IHG22','IHG23','IHG24','ACD1','ACD2','ACD3','ACD4','ACD5','ACD6','PCD1','PCD2','PCD3','PCD4','PCD5','PCD6','SFS1','SFS2','SFS3','SFS4','IFS2','IFS3','IFS4','OFS3','OFS4','ATS5','ATS6','ATS7','ATS8','PTS5','PTS6','PTS7','PTS8'};
labels(9).focus    = {'LFG3','LFG4','LFG5','LFG6','LFG7','LFG10','LFG11','LFG12','LFG13','LFG14','LFG15','LFG19','LFG20','LFG21','LFG22','LFG23','LFG24','LFG27','LFG28','LFG29','LFG30','LFG31','LFG32','LFG35','LFG36','LFG37','LFG38','LFG39','LFG40','LFG44','LFG45','LFG46','LFG47','LFG48','LFG56','LFG57','IHG1','IHG2','IHG3','IHG4','IHG5','IHG6','IHG9','IHG10','IHG11','IHG12','IHG13','IHG14'};
labels(9).channel  = [2 4:14 16:31 33:106 108:117];

labels(10).subject = 'PY11N006';
labels(10).values  = {'LAF2','LAF3','LAF4','LOF1','LOF2','LOF3','LOF4','LOF5','LOF6','LOF7','LOF8','LFP1','LFP2','LFP3','LFP4','LAT1','LAT3','LAT4','LAT5','LAT6','LMT1','LMT2','LMT3','LMT4','LMT5','LMT6','LPT1','LPT2','LPT3','LPT4','LPT5','LPT6','RAF2','RAF3','RAF4','ROF1','ROF2','ROF3','ROF4','ROF5','ROF6','ROF7','ROF8','RFP1','RFP2','RFP3','RFP4','RAT1','RAT2','RAT3','RAT4','RAT5','RAT6','RMT1','RMT2','RMT3','RMT4','RMT5','RMT6','RPT1','RPT2','RPT3','RPT4','RPT5','RPT6'};
labels(10).focus   = {'RAT1','RAT2','RAT3','RAT4','RAT5','RAT6','RMT1','RMT2','RMT3','RMT4','RMT5','RMT6'};
labels(10).channel = [1:16 18:66];

labels(11).subject = 'PY12N005';
labels(11).values  = {'LLT4','LLT5','LLT6','LLT7','LLT8','ABT1','ABT2','ABT3','ABT4','ABT5','ABT6','PBT1','PBT2','PBT3','PBT4','PBT5','PBT6','ADL1','ADL2','ADL3','ADL4','ADL5','ADL6','ADL7','ADL8','PDL1','PDL2','PDL3','PDL4','PDL5','PDL6','PDL7','PDL8','AHD1','AHD2','AHD3','AHD4','AHD5','AHD6','AHD7','AHD8','LPG1','LPG2','LPG3','LPG4','LPG5','LPG6','LPG7','LPG8','LPG9','LPG10','LPG11','LPG12','LPG13','LPG14','LPG15','LPG16','LPG17','LPG18','LPG19','LPG20','LPG21','LPG22','LPG23','LPG24','LPG25','LPG26','LPG27','LPG28','LPG29','LPG30','LPG31','LPG32','LPG33','LPG34','LPG35','LPG36','LPG37','LPG38','LPG39','LPG40','LPG41','LPG42','LPG43','LPG44','LPG45','LPG46','LPG47','LPG48','LPG49','LPG50','LPG51','LPG52','LPG53','LPG54','LPG55','LPG56','LPG57','LPG58','LPG59','LPG60','LPG61','LPG62','LPG63','LPG64'};
labels(11).focus   = {'LPG25','LPG26','LPG27','LPG28','LPG29','LPG33','LPG34','LPG35','LPG36','LPG37','LPG41','LPG42','LPG43','LPG44','LPG45','LPG46','LPG49','LPG50','LPG51','LPG52','LPG53','LPG54','LPG58','LPG59','LPG60','LPG61','LPG62'};
labels(11).channel = [1:4 7:107];

labels(12).subject = 'PY12N008';
labels(12).values  = {'AFS1','AFS2','AFS3','AFS4','AFS5','AFS6','AFS7','AFS8','OFS1','OFS2','OFS3','OFS4','OFS5','OFS6','OFS7','OFS8','ATS1','ATS2','ATS3','ATS4','BTS1','BTS2','BTS3','BTS4','BTS5','BTS6','FPG1','FPG2','FPG3','FPG4','FPG5','FPG6','FPG9','FPG10','FPG11','FPG12','FPG13','FPG14','FPG15','FPG17','FPG18','FPG19','FPG20','FPG21','FPG22','FPG23','FPG24','FPG25','FPG26','FPG27','FPG28','FPG29','FPG30','FPG31','FPG32','FPG33','FPG34','FPG35','FPG36','FPG37','FPG38','FPG39','FPG40','FPG41','FPG42','FPG43','FPG44','FPG45','FPG46','FPG47','FPG48','FPG49','FPG50','FPG51','FPG52','FPG53','FPG54','FPG55','FPG56','FPG57','FPG58','FPG59','FPG60','FPG61','FPG62','FPG63','FPG64'};
labels(12).focus   = {'FPG34','FPG35','FPG41','FPG42','FPG43','FPG44','FPG49','FPG50','FPG51','FPG52','FPG57','FPG58','FPG59','FPG60'};
labels(12).channel = [1:4 7:89];

labels(13).subject = 'PY12N010';
labels(13).values  = {'AFS2','AFS3','AFS4','AFS5','AFS1','MFS1','AFS6','MFS3','MFS4','MFS5','MFS6','MFS7','MFS8','ATS1','ATS2','ATS3','ATS4','ATS5','ATS6','ABS1','ABS2','ABS3','ABS4','ABS5','ABS6','PBS1','PBS2','PBS3','PBS4','PBS5','PBS6','LAD1','LAD2','LAD3','LAD4','LAD5','LAD6','LAD7','LHD1','LHD2','LHD3','LHD4','LHD6','LHD7','LHD8','FTG1','FTG2','FTG3','FTG4','FTG5','FTG6','FTG7','FTG8','FTG9','FTG10','FTG11','FTG12','FTG13','FTG14','FTG15','FTG16','FTG17','FTG18','FTG19','FTG20','FTG21','FTG22','FTG23','FTG24','FTG25','FTG26','FTG27','FTG28','FTG29','FTG30','FTG31','FTG32','FTG33','FTG34','FTG35','FTG36','FTG37','FTG38','FTG39','FTG40','FTG41','FTG42','FTG43','FTG44','FTG45','FTG46','FTG47','FTG48','FTG49','FTG50','FTG51','FTG52','FTG53','FTG54','FTG55','FTG56','FTG57','FTG58','FTG59','FTG60','FTG61','FTG62','FTG63','FTG64'};
labels(13).focus   = {'ABS1','ABS2','ABS3','ABS4','ABS5','ABS6','PBS1','PBS2','PBS3','PBS4','PBS5','PBS6','FTG49','FTG50','FTG51','FTG52','FTG53','FTG54','FTG57','FTG58','FTG59','FTG60','FTG61','FTG62'};
labels(13).channel = [1:38 40:43 45:111];

labels(14).subject = 'PY12N012';
labels(14).values  = {'LIF1','LIF2','LIF3','LIF4','LIF5','LIF6','LOF2','LOF3','LOF4','LOF5','LOF6','LAT1','LAT2','LAT3','LAT4','ABT1','ABT2','ABT3','ABT4','MBT1','MBT2','MBT3','MBT4','PBT1','PBT2','PBT3','PBT4','PBT5','PBT6','LAD1','LAD2','LAD3','LAD4','LAD5','LAD6','LAD7','LAD8','LHD1','LHD2','LHD3','LHD4','LHD5','LHD6','LHD7','LHD8','LFT1','LFT2','LFT3','LFT4','LFT5','LFT6','LFT7','LFT8','LFT9','LFT10','LFT11','LFT12','LFT13','LFT14','LFT15','LFT16','LFT17','LFT18','LFT19','LFT20','LFT21','LFT22','LFT23','LFT24','LFT25','LFT26','LFT27','LFT28','LFT29','LFT30','LFT31','LFT32','LFT33','LFT34','LFT35','LFT36','LFT37','LFT38','LFT39','LFT40','LFT41','LFT42','LFT43','LFT44','LFT45','LFT46','LFT47','LFT48','LFT49','LFT50','LFT51','LFT52','LFT53','LFT54','LFT55','LFT56','LFT57','LFT58','LFT59','LFT60','LFT61','LFT62','LFT63','LFT64'};
labels(14).focus   = {'LFT1','LFT2','LFT3','LFT4','LFT5','LFT9','LFT10','LFT11','LFT12','LFT13','LFT14','LFT17','LFT18','LFT19','LFT20','LFT21','LFT22','LFT25','LFT26','LFT27','LFT28','LFT29','LFT33','LFT34','LFT35','LFT36'};
labels(14).channel = 1:109;

labels(15).subject = 'PY13N001';
labels(15).values  = {'OCIS1','OCIS2','OCIS3','OCIS4','OCIS17','OCIS32','OCIS5','OCIS6','OCIS7','OCIS8','OCIS9','OCIS10','OCIS11','OCIS12','OCIS13','OCIS14','OCIS15','OCIS16','OCIS18','OCIS19','OCIS20','OCIS21','OCIS22','OCIS23','OCIS24','OCIS25','OCIS26','OCIS27','OCIS28','OCIS29','OCIS30','OCIS31','LMOS1','LMOS2','LMOS3','LMOS4','LMOS5','LMOS6','LMOS7','LMOS8','LLOS1','LLOS2','LLOS3','LLOS4','LLOS5','LLOS6','LAOS1','LAOS2','LAOS3','LAOS4','LAOS5','LAOS6','LAOS7','LAOS8','LLTS1','LLTS2','LLTS3','LLTS4','LLTS5','LLTS6','LLFS1','LLFS2','LLFS3','LLFS4','LOPG1','LOPG2','LOPG3','LOPG4','LOPG5','LOPG6','LOPG7','LOPG8','LOPG9','LOPG10','LOPG11','LOPG12','LOPG13','LOPG14','LOPG15','LOPG16','LOPG17','LOPG18','LOPG19','LOPG20','LOPG21','LOPG22','LOPG23','LOPG24','LOPG25','LOPG26','LOPG27','LOPG28','LOPG29','LOPG30','LOPG31','LOPG32','LOPG33','LOPG34','LOPG35','LOPG36','LOPG37','LOPG38','LOPG39','LOPG40','LOPG41','LOPG42','LOPG43','LOPG44','LOPG45','LOPG46','LOPG47','LOPG48','LOPG49','LOPG50','LOPG51','LOPG52','LOPG53','LOPG54','LOPG55','LOPG56','LOPG57','LOPG58','LOPG59','LOPG60','LOPG61','LOPG62','LOPG63','LOPG64'};
labels(15).focus   = {'LMOS8','LLOS6','LOPG9','LOPG10','LOPG17','LOPG18','LOPG19','LOPG25','LOPG26','LOPG27','LOPG33','LOPG34','LOPG35','LOPG41','LOPG42','LOPG43','LOPG49','LOPG50','LOPG51','LOPG57','LOPG58'};
labels(15).channel = 1:128;

labels(16).subject = 'PY13N003';
labels(16).values  = {'RTO1','RTO2','RTO3','RTO4','RTO5','RTO6','RTO7','RTO8','RTO9','RTO10','RTO11','RTO12','RTO13','RTO14','RTO15','RTO16','RTO17','RTO18','RTO19','RTO20','RAT1','RAT2','RAT3','RAT4','RMT1','RMT2','RMT3','RMT4','RMT5','RMT6','RPT1','RPT2','RPT3','RPT4','RPT5','RPT6','RPT7','RPT8','ROF1','ROF2','ROF3','ROF4','ROF5','ROF6','RIF1','RIF2','RIF3','RIF4','RIF5','RIF6','RFT1','RFT2','RFT3','RFT4','RFT5','RFT6','RFT7','RFT8','RFT9','RFT10','RFT11','RFT12','RFT13','RFT14','RFT15','RFT16','RFT17','RFT18','RFT19','RFT20','RFT21','RFT22','RFT23','RFT24','RFT25','RFT26','RFT27','RFT28','RFT29','RFT30','RFT31','RFT32','RFT33','RFT34','RFT35','RFT36','RFT37','RFT38','RFT39','RFT40','RFT41','RFT42','RFT43','RFT44','RFT45','RFT46','RFT47','RFT48','RFT49','RFT50','RFT51','RFT52','RFT53','RFT54','RFT55','RFT56','RFT57','RFT58','RFT59','RFT60','RFT61','RFT62','RFT63','RFT64','ALD1','ALD2','ALD3','ALD4','ALD5','ALD6','ALD7','ALD8','PLD1','PLD2','PLD3','PLD4','PLD5','PLD6','PLD7','PLD8','RHD1','RHD2','RHD3','RHD4','RHD5','RHD6','RHD7','RHD8'};
labels(16).focus   = {'RTO10','RTO15','RTO20','RMT1','RMT2','RMT3','RMT4','RMT5','RMT6','RPT1','RPT2','RPT3','RPT4','RPT5','RPT6','RPT7','RPT8','RFT13','RFT14','RFT15','RFT21','RFT22','RFT23','RFT24','RFT31','RFT32','RFT40','ALD1','ALD2','ALD3','ALD4','ALD5','ALD6','ALD7','ALD8','PLD1','PLD2','PLD3','PLD4','PLD5','PLD6','PLD7','PLD8'};
labels(16).channel = [1:4 7:140];

labels(17).subject = 'PY13N004';
labels(17).values  = {'TPO1','TPO2','TPO3','TPO4','TPO7','TPO8','TPO9','TPO10','TPO11','TPO12','TPO13','TPO14','TPO15','TPO16','TPO17','TPO18','TPO19','TPO20','TPO21','TPO22','TPO23','TPO24','TPO25','TPO26','TPO27','TPO28','TPO29','TPO30','TPO31','TPO32','TPO33','TPO34','TPO35','TPO36','TPO37','TPO38','TPO39','TPO40','TPO41','TPO42','TPO43','TPO44','TPO45','TPO46','TPO47','TPO48','TPO49','TPO50','TPO51','TPO52','TPO53','TPO54','TPO55','TPO56','TPO57','TPO58','TPO59','TPO60','TPO61','TPO62','TPO63','TPO64','ILD1','ILD2','ILD3','ILD4','ILD5','ILD6','ILD7','ILD8','RTL1','RTL2','RTL3','RTL4','RTL5','RTL6','RTL7','RTL8','PBT1','PBT2','PBT3','PBT4','PBT5','PBT6','ABO1','ABO2','ABO3','ABO4','ABO5','ABO6','PBO1','PBO2','PBO3','PBO4','SLD1','SLD2','SLD3','SLD4','SLD5','SLD6','SLD7','SLD8','ABT1','ABT2','ABT3','ABT4','ABT5','ABT6','ABT7','ABT8','RIH1','RIH2','RIH3','RIH4','RIH5','RIH6','RIH7','RIH8','RIH9','RIH10','RIH11','RIH12','RIH13','RIH14','RIH15','RIH16','LIH17','LIH18','LIH19','LIH20','LIH21','LIH22','LIH23','LIH24','LIH25','LIH26','LIH27','LIH28','LIH29','LIH30','LIH31','LIH32'};
labels(17).focus   = {'TPO33','TPO34','TPO35','TPO41','TPO42','TPO43','TPO49','TPO50','TPO51','TPO57','TPO58','TPO59','ABO1','ABO2','ABO3','ABO4','ABO5','ABO6','PBO1','PBO2','PBO3','PBO4','LIH17','LIH18','LIH19','LIH20','LIH21','LIH22','LIH23','LIH24','LIH25','LIH26','LIH27','LIH28','LIH29','LIH30','LIH31','LIH32'};
labels(17).channel = [1:4 7:144];

labels(18).subject = 'PY13N010';
labels(18).values  = {'IFS1','IFS2','IFS3','IFS4','PBT1','PBT2','PBT3','PBT4','PBT5','PBT6','SFS1','SFS2','SFS3','SFS4','SFS5','SFS6','IPS1','IPS2','IPS3','IPS4','ABT1','ABT2','ABT3','ABT4','MBT1','MBT2','MBT3','MBT4','MFS1','MFS2','MFS3','MFS4','MFS5','MFS6','LAD1','LAD2','LAD3','LAD4','LAD5','LAD6','LAD7','LAD8','LHD1','LHD2','LHD3','LHD4','LHD5','LHD6','LHD7','LHD8','TPG1','TPG2','TPG3','TPG4','TPG5','TPG6','TPG7','TPG8','TPG9','TPG10','TPG11','TPG12','TPG13','TPG14','TPG15','TPG16','TPG17','TPG18','TPG19','TPG20','TPG21','TPG22','TPG23','TPG24','TPG25','TPG26','TPG27','TPG28','TPG29','TPG30','TPG31','TPG32','FPG33','FPG34','FPG35','FPG36','FPG37','FPG38','FPG39','FPG40','FPG41','FPG42','FPG43','FPG44','FPG45','FPG46','FPG47','FPG48','FPG49','FPG50','FPG51','FPG52','FPG53','FPG54','FPG55','FPG56','FPG57','FPG58','FPG59','FPG60','FPG61','FPG62','FPG63','FPG64'};
labels(18).focus   = {'ABT1','ABT2','ABT3','ABT4','LAD1','LAD2','LAD3','LAD4','LAD5','LAD6','LAD7','LAD8','LHD1','LHD2','LHD3','LHD4','LHD5','LHD6','LHD7','LHD8','TPG1','TPG2','TPG3','TPG4','TPG5','TPG9','TPG10','TPG11','TPG12','TPG13','TPG17','TPG18','TPG19','TPG20'};
labels(18).channel = [1:4 7:116];

labels(19).subject = 'PY13N011';
labels(19).values  = {'RBT1','RBT2','RBT3','RBT4','RBT5','RBT6','RBT7','RBT8','RBO1','RBO2','RBO3','RBO4','RBO5','RBO6','RTO1','RTO2','RTO3','RTO4','RTO5','RTO6','RTO7','RTO8','RTO9','RTO10','RTO11','RTO12','RTO13','RTO14','RTO15','RTO16','RTO17','RTO18','RTO19','RTO20','RTG1','RTG2','RTG3','RTG4','RTG5','RTG6','RTG7','RTG8','RTG9','RTG10','RTG11','RTG12','RTG13','RTG14','RTG15','RTG16','RTG17','RTG18','RTG19','RTG20','RTG21','RTG22','RTG23','RTG24','RTG25','RTG26','RTG27','RTG28','RTG29','RTG30','RTG31','RTG32','RTG33','RTG34','RTG35','RTG36','RTG37','RTG38','RTG39','RTG40','RTG41','RTG42','RTG43','RTG44','RTG45','RTG46','RTG47','RTG48'};
labels(19).focus   = {'RBO5','RBO6','RTO4','RTO5','RTO9','RTO10','RTO14','RTO15','RTO19','RTO20','RTG4','RTG5','RTG6','RTG7','RTG8','RTG12','RTG13','RTG14','RTG15','RTG16','RTG20','RTG21','RTG22','RTG23','RTG24','RTG28','RTG29','RTG30','RTG31','RTG32','RTG36','RTG37','RTG38','RTG39','RTG40','RTG44','RTG45','RTG46','RTG47','RTG48'};
labels(19).channel = [1:4 7:84];

labels(20).subject = 'PY14N004';
labels(20).values  = {'RAT1','RAT2','RAT3','RAT4','RIT1','RIT2','RIT3','RIT4','RPT1','RPT2','RPT3','RPT4','RPT5','RPT6','RST1','RST3','RST4','RHD1','RHD2','RHD3','RHD4','RHD5','RHD6','RAD1','RAD2','RAD3','RAD4','RAD5','RAD6','RAT5','RAT6','LAT1','LAT2','LAT3','LAT4','LAT5','LAT6','LIT1','LIT2','LIT3','LIT4','LPT1','LPT2','LPT3','LPT4','LPT5','LPT6','LST1','LST2','LST3','LST4','LHD1','LHD2','LHD3','LHD4','LHD5','LHD6','LHD7','LHD8','LAD1','LAD2','LAD3','LAD4','LAD5','LAD6','LAD7','LAD8'};
labels(20).focus   = {'LAT1','LAT2','LAT3','LAT4','LAT5','LAT6','LIT1','LIT2','LIT3','LIT4','LHD1','LHD2','LHD3','LHD4','LHD5','LHD6','LHD7','LHD8','LAD1','LAD2','LAD3','LAD4','LAD5','LAD6','LAD7','LAD8'};
labels(20).channel = [1:15 17:24 27:70];

labels(21).subject = 'PY14N005';
labels(21).values  = {'PFD1','PFD2','PFD3','PFD4','PFD5','PFD6','PFD7','PFD8','AFD1','AFD2','AFD3','AFD4','AFD5','AFD6','AFD7','AFD8','SFD1','SFD2','SFD3','SFD4','SFD5','SFD6','SFD7','SFD8','FTG1','FTG2','FTG3','FTG4','FTG5','FTG6','FTG7','FTG8','FTG9','FTG10','FTG11','FTG12','FTG13','FTG14','FTG15','FTG16','FTG17','FTG18','FTG19','FTG20','FTG21','FTG22','FTG23','FTG24','FTG25','FTG26','FTG27','FTG28','FTG29','FTG30','FTG31','FTG32','FTG33','FTG34','FTG35','FTG36','FTG37','FTG38','FTG39','FTG40','FTG41','FTG42','FTG43','FTG44','FTG45','FTG46','FTG47','FTG48','FTG49','FTG50','FTG51','FTG52','FTG53','FTG54','FTG55','FTG56','FTG57','FTG58','FTG59','FTG60','FTG61','FTG62','FTG63','FTG64'};
labels(21).focus   = {'PFD1','PFD2','PFD3','PFD4','PFD5','PFD6','PFD7','PFD8','AFD1','AFD2','AFD3','AFD4','AFD5','AFD6','AFD7','AFD8','SFD1','SFD2','SFD3','SFD4','SFD5','SFD6','SFD7','SFD8','FTG20','FTG21','FTG22','FTG23','FTG28','FTG29','FTG30','FTG31','FTG36','FTG37','FTG38','FTG39','FTG44','FTG45','FTG46','FTG47'};
labels(21).channel = 1:88;

labels(22).subject = 'PY15N004';
labels(22).values  = {'IFS3','IFS4','IFS5','IFS6','IFS1','IFS2','IFS7','IFS8','ATS1','ATS2','ATS5','ATS6','MTS2','MTS3','MTS4','MTS5','MTS6','PTS1','PTS2','PTS3','ATS3','ATS4','PTS4','PTS5','PTS6','RAD1','RAD2','RAD3','RAD4','RAD5','RAD6','RAD7','RAD8','AHD1','AHD2','AHD3','AHD4','AHD5','AHD6','AHD7','AHD8','PHD1','PHD2','PHD3','PHD4','PHD5','PHD6','PHD7','PHD8','FTG1','FTG2','FTG3','FTG4','FTG5','FTG6','FTG7','FTG8','FTG9','FTG10','FTG11','FTG12','FTG13','FTG14','FTG15','FTG16','FTG17','FTG18','FTG19','FTG20','FTG21','FTG22','FTG23','FTG24','FTG25','FTG26','FTG27','FTG28','FTG29','FTG30','FTG31','FTG32','FTG33','FTG34','FTG35','FTG36','FTG37','FTG38','FTG39','FTG40','FTG41','FTG42','FTG43','FTG44','FTG45','FTG46','FTG47','FTG48','FTG49','FTG50','FTG51','FTG52','FTG53','FTG54','FTG55','FTG56','FTG57','FTG58','FTG59','FTG60','FTG61','FTG62','FTG63','FTG64'};
labels(22).focus   = {};
labels(22).channel = [1:12 14:19 21:37 42 43 46:57 60:123];


% store the struct array in a .mat file
save('infolabels.mat','labels');

fprintf('The End\n');