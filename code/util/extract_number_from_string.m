function s_idx = extract_number_from_string(s_name)
s_name(strfind(s_name, '_')) = [];
s_idx = s_name;
deleteMe = isletter(s_name);
s_idx(deleteMe) = [];

end