% This is a pipeline to process ISS example images all way from the
% beginning, i.e. raw acquisition files.
% The provided files are in Nikon .nd2 format, which supports bio-formats
% Cartana, Xiaoyan, 2020



% Initiate
ctn = cartana;
ctn.OutputDirectory = 'W:\Analysis\Christin Mueller\ISS\exp64\236KS\236KS_AP1_5cycles';

% How many rounds of ISS cycles
ctn.nRounds = 5;

% EDIT: =1 if background channels present else 0
ctn.BackgroundChannels = 1;
ctn.extraChannel = 1;

% Whether library prep image is provided and required to be used as reference
% If no library preparation image was taken, the value is 0, and ISS cycle 1
% images will be used as reference instead
ctn.isLibraryPrep = 0;

% The channel number of anchor stain in library prep image
% If library prep image anchor and DAPI images are provided, this image
% will be used as reference ISS spots
ctn.AnchorChannel = 2;

% DAPI channel, one for each ISS cycle, Lirary Prep is the first
ctn.DAPIChannel = [1 1 1 1 1 1 1]; 

% ISS channels
% One row for one cycle of ISS, 4 channels in each cycle
% Cycle order: AlexFluor750, Cy3, Cy5, AlexFluor488 (same in all cycles)
ctn.ISSChannels = repmat([5 3 4 2], 5 + ctn.BackgroundChannels + ctn.extraChannel, 1);    


%% Pre-align stitched images 
% To compensate for inaccuracy in microscope stage repositioning 
% EDIT: add background channels as last Round(line) and set ctn.BackgroundChannels to 1


ctn.StitchedImages = {...
'W:\Analysis\Christin Mueller\ISS\exp64\236KS\236KS_AP1_DAPI.tif', 'W:\Analysis\Christin Mueller\ISS\exp64\236KS\236KS_AP1_FITC.tif', 'W:\Analysis\Christin Mueller\ISS\exp64\236KS\236KS_AP1_Cy3.tif', 'W:\Analysis\Christin Mueller\ISS\exp64\236KS\236KS_AP1_Cy5.tif', 'W:\Analysis\Christin Mueller\ISS\exp64\236KS\236KS_AP1_Cy7.tif';
'W:\Analysis\Christin Mueller\ISS\exp64\236KS\236KS_AP2_DAPI.tif', 'W:\Analysis\Christin Mueller\ISS\exp64\236KS\236KS_AP2_FITC.tif', 'W:\Analysis\Christin Mueller\ISS\exp64\236KS\236KS_AP2_Cy3.tif', 'W:\Analysis\Christin Mueller\ISS\exp64\236KS\236KS_AP2_Cy5.tif', 'W:\Analysis\Christin Mueller\ISS\exp64\236KS\236KS_AP2_Cy7.tif';
'W:\Analysis\Christin Mueller\ISS\exp64\236KS\236KS_AP3_DAPI.tif', 'W:\Analysis\Christin Mueller\ISS\exp64\236KS\236KS_AP3_FITC.tif', 'W:\Analysis\Christin Mueller\ISS\exp64\236KS\236KS_AP3_Cy3.tif', 'W:\Analysis\Christin Mueller\ISS\exp64\236KS\236KS_AP3_Cy5.tif', 'W:\Analysis\Christin Mueller\ISS\exp64\236KS\236KS_AP3_Cy7.tif';
'W:\Analysis\Christin Mueller\ISS\exp64\236KS\236KS_AP4_DAPI.tif', 'W:\Analysis\Christin Mueller\ISS\exp64\236KS\236KS_AP4_FITC.tif', 'W:\Analysis\Christin Mueller\ISS\exp64\236KS\236KS_AP4_Cy3.tif', 'W:\Analysis\Christin Mueller\ISS\exp64\236KS\236KS_AP4_Cy5.tif', 'W:\Analysis\Christin Mueller\ISS\exp64\236KS\236KS_AP4_Cy7.tif';
'W:\Analysis\Christin Mueller\ISS\exp64\236KS\236KS_AP5_DAPI.tif', 'W:\Analysis\Christin Mueller\ISS\exp64\236KS\236KS_AP5_FITC.tif', 'W:\Analysis\Christin Mueller\ISS\exp64\236KS\236KS_AP5_Cy3.tif', 'W:\Analysis\Christin Mueller\ISS\exp64\236KS\236KS_AP5_Cy5.tif', 'W:\Analysis\Christin Mueller\ISS\exp64\236KS\236KS_AP5_Cy7.tif';
'W:\Analysis\Christin Mueller\ISS\exp64\236KS\236KS_BG_DAPI.tif', 'W:\Analysis\Christin Mueller\ISS\exp64\236KS\236KS_BG_FITC.tif', 'W:\Analysis\Christin Mueller\ISS\exp64\236KS\236KS_BG_Cy3.tif', 'W:\Analysis\Christin Mueller\ISS\exp64\236KS\236KS_BG_Cy5.tif', 'W:\Analysis\Christin Mueller\ISS\exp64\236KS\236KS_BG_Cy7.tif';
'W:\Analysis\Christin Mueller\ISS\exp64\236KS\236KS_AP1_DAPI.tif', 'W:\Analysis\Christin Mueller\ISS\exp64\236KS\236KS_AP1_FITC.tif', 'W:\Analysis\Christin Mueller\ISS\exp64\236KS\236KS_AP1_Cy3.tif', 'W:\Analysis\Christin Mueller\ISS\exp64\236KS\236KS_AP1_Cy5.tif', 'W:\Analysis\Christin Mueller\ISS\exp64\236KS\236KS_AP1_Cy7.tif';

};


                                                                                                                                                                                                                                                                                                                                          
if ctn.BackgroundChannels
    ctn.nRounds = ctn.nRounds + 1;
end

if ctn.extraChannel
    ctn.nRounds = ctn.nRounds + 1;
end

ctn = prealign_sections(ctn);

%% Tile images into a size that is easy to handle
% These can also be input files for starfish pipeline
ctn = tile_images(ctn);

if ctn.BackgroundChannels
    ctn.nRounds = ctn.nRounds - 1;
end

if ctn.extraChannel
    ctn.nRounds = ctn.nRounds - 1;
end

%% Prepare matrix to specify file 

if ctn.BackgroundChannels
    ctn = specify_input_edit_extraChannel(ctn);
else
    ctn = specify_input(ctn);
end
save(fullfile(ctn.OutputDirectory, 'FilesReady.mat'), 'ctn');

%% Use Cellprofiler pipeline to 
% Run cellprofiler using the InputForCP.csv in the first module

% decodeing
ctn.TaglistFile = '‪W:\A=nalysis\Christin Mueller\codebooks\codebook_al_5cycles_1_2_3_4_5.csv'
ctn.SaveIntermediate=1;
ctn = decode_cellprofiler(ctn);

% save MATLAB variable file, but may take a long time
save(fullfile(ctn.OutputDirectory, 'CP.mat'), 'ctn', '-v7.3');

% % or use MATLAB
% ctn.MinSignalIntensity = 400; % minimal intensity level in the original image to be considered above background
% ctn = decode_simplified_inmatlab(ctn);

%% plot againplot_reads(ctn)
