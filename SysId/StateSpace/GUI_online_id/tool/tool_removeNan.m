function cleanedRows = tool_removeNan(rows)

cleanedRows = rows;
%if there is a Nan in a row, delte that row from rows, at the end print how
%many rows have been deleted
deletedRows = 0;

%take every columns
col_1 = rows{1};
col_2 = rows{2};
col_3 = rows{3};
col_4 = rows{4};

%index where nan values are
indexNan = [];

%every col mus have the same length
N = length(col_1);

%see if in row i there is at least one Nan value
for i = 1 : N
   if(isnan(col_1(i)) || isnan(col_2(i)) || isnan(col_3(i)) || isnan(col_4(i)))
       %there is at least one Nan value, update indexNan and deletedRows
       indexNan = [indexNan, i];
       deletedRows = deletedRows + 1;
   end
end

if(deletedRows > 0)
    %delte rows
    cleanedRows{1}(indexNan) = [];
    cleanedRows{2}(indexNan) = [];
    cleanedRows{3}(indexNan) = [];
    cleanedRows{4}(indexNan) = [];
    %print how many rows have been deleted
    display(['Deleted ' num2str(deletedRows) ' row(s) which had a Nan value']);
end

end

