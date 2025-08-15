

function dataset = extract_dataset(animal_name, area_name, options)

    
%     switch animal_name
%         case 'chip'
%             visual_delay = 200;
%             photo_channel = 7;
%         case 'hans'
%             visual_delay = 250;
%             photo_channel = 7;
%         otherwise
%             error('Wrong animal');
%     end

    dataset = localization_extract_data(animal_name, area_name);
    dataset = localization_eyeposition(dataset, 'eye', 1);

end

    
