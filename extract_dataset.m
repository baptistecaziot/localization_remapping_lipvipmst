

function dataset = extract_dataset(animal_name, area_name, settings)


    arguments
        animal_name {mustBeText}
        area_name {mustBeText}
        settings {mustBeA(settings,"struct")}
    end

    dataset = localization_extract_data(animal_name, area_name, settings);
    dataset = localization_eyeposition(dataset, settings);

end

    
