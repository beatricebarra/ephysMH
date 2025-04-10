function replaced_cell_array = replace_elements_cellarray(input_cell_array, condition, cond_value, replace_value)
    replaced_array = input_cell_array;
    eval( ['index = find(input_cell_array', condition,  num2str(cond_value), ')']);
    replaced_array(index) = replace_value;
end