function RevInd = FindRevs(start,stop,g10)
RevInd = zeros(length(start),1);
for kk=1:length(start)
    if sign(g10(start(kk))) ~= sign(g10(stop(kk)))
        RevInd(kk) = 1;
    end
end
