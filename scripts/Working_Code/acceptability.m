function p = acceptability(V, F)
%ACCEPTABILITY Evaluates percent of faces with aspect ratio under 3

AR = zeros(length(F),1);
for i=1:length(F)
    AR(i) = aspect_ratio(V(F(i,:),:));
end

p = sum(AR < 3) / length(AR);

    function ar = aspect_ratio(triangle)
        a = norm(triangle(2,:) - triangle(1,:));
        b = norm(triangle(3,:) - triangle(2,:));
        c = norm(triangle(1,:) - triangle(3,:));
        
        s = (a + b + c) / 2;
        
        ar = (a*b*c) / (8*(s-a)*(s-b)*(s-c));
    end

end

