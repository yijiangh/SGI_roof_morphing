function c = complexity(V, F)
%COMPLEXITY returns "complexity" of a shape--perimeter/area

O = outline(F);
perimeter = sum(normrow(V(O(:,2),:) - V(O(:,1),:)));

area = sum(doublearea(V, F)) / 2;

c = perimeter / area;

end

