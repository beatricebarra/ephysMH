function [maxidx]= find_true_peak(FR, xmax)
if isempty(xmax.loc)
    maxidx = NaN; 
    %disp('case0')
else  
%     if xmax.loc(1)==1 && length(xmax.loc) ==1% if it's the first element and the vector has only one element discard
%         maxidx = NaN; 
%         disp('case1')
    if FR(xmax.loc(1))< FR(xmax.loc(1)-1) % negative derivative, means this is not a peak
        if length(xmax.loc)>1 % are there more elements, if yes...
            maxidx = xmax.loc(2); % Take second peak
        else % if not... 
            maxidx = NaN; % NaN
        end
    elseif xmax.loc(1)==length(FR) % if it's the last element discard
         maxidx = NaN; 
         %disp('case2')
%     elseif xmax.loc(1)==1 && length(xmax.loc)>1 % if it's the first element and the vector has more then one element, take second peak
%         if xmax.loc(2)~=length(FR)
%             maxidx = xmax.loc(2); 
%         else
%             maxidx = NaN; % if the second element is also the last one, discard
%             disp('case3')
%         end
    else
        maxidx = xmax.loc(1); % good case if which the first peak is also ok
    end
    
end