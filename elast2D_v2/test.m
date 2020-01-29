a = [1; 2; 3; 4 ;5 ;6; 7; 8; 9];
b = a(randperm(length(a),3));
c = [];
for i = 1:length(a);
    if a(i) ~= b;
        c = [ c ;
             a(i)];
    end
end
b
c