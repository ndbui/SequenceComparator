"""
    append_dict!(source::Dict, new_data::Dict)

Returns the source dictionary with key, value pairs from the new_data dictionary appended to it
"""
function append_dict!(source::Dict, new_data::Dict)
    for (key, value) in new_data
        if !(key in keys(source))
            source[key] = value
        end
    end
    return source
end