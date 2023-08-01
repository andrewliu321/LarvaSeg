%%To measure length, width, area, volume of larvas

folder = 'C:\Users\Andrew\Desktop\LarvaSeg\';
% this is where you're image files exist


% File Names:
genotype = 'Example Larva ';
% usually the larva has a genotype such as WT or mutants and might have age
% followed by sample number


%chose which file to load:                   
for position =1
    %if running through a group of samples, do position = 1:10
    

    % run set up if first time
    if position == 1

    cd 'C:\Users\Andrew\Desktop\LarvaSeg\'
    set(0,'DefaultFigureWindowStyle','docked')
    end

    %%
    im = imread([folder,genotype,num2str(position),'.tif']);

    f1 = figure(1)
    imshow(im)

    %% hand segment
    draw_larva
    %%
    %record_mid_index = new points identifying the midline (pixel list)
    %larva_record = record of segmented larva


    % measure pairwise distance lengths
    clear larva_length_per_pix
    clear xx
    for i = 1:length(record_mid_index{position})-1
        xx = record_mid_index{position}(i:i+1,:);
        larva_length_per_pix(i)=pdist(xx,'euclidean');
    end

    larva_length(position) = sum(larva_length_per_pix)
    %records the length of the midline drawn
    
    larva_avg_width(position) = mean(nonzeros(d))
    %averages the width measured in orthogonal lines
    
    larva_volume(position) = larva_length(position)*pi*(larva_avg_width(position)/2)^2
    %estimates the volume using the average width/2 as radius, length and
    %cone formula
    
    larva_area(position)=bwarea(larva2)
    %calculate the area of the segmented larva

    %all numbers are in pixels, need to convert to um separately.
end




%% Save data

save_file_name = ['Measured_Larva',genotype,'.mat'];

save(save_file_name, 'record_mid_index','larva_record',...
    'larva_length','larva_avg_width','larva_volume','larva_area')



transpose(larva_length)

transpose(larva_avg_width)

transpose(larva_volume)

transpose(larva_area)









