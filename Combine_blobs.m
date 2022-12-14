%% 
images ={...
    'W:\Analysis\Christin Mueller\ISS\exp64\233KS\233KS_AP1\files 1\blobs.csv',
    'W:\Analysis\Christin Mueller\ISS\exp64\233KS\233KS_AP1\files 2\blobs.csv',
    'W:\Analysis\Christin Mueller\ISS\exp64\233KS\233KS_AP1\files 3\blobs.csv',
    'W:\Analysis\Christin Mueller\ISS\exp64\233KS\233KS_AP1\files 4\blobs.csv',
    'W:\Analysis\Christin Mueller\ISS\exp64\233KS\233KS_AP1\files 5\blobs.csv',
    'W:\Analysis\Christin Mueller\ISS\exp64\233KS\233KS_AP1\files 6\blobs.csv',

  };

imageindex = {1 170 339 508 677 846};

data= readtable(images{1});
for i =1:length(images)-1
    
    dataadd= readtable(images{i+1});
    dataadd{:,1}=dataadd{:,1}+imageindex{i+1}-1;
    data = vertcat(data,dataadd);
end   
writetable(data,'W:\Analysis\Christin Mueller\ISS\exp64\233KS\blobs.csv');
%% 