function replaced_cell_array = fill_empty_elements_cellarray(input_cell_array, replace_value)
    replaced_array = input_cell_array;
    index = isempty(input_cell_array); 
    replaced_array(index) = replace_value;
end