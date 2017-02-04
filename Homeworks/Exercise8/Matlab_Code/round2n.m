function [rounded] = round2n(in,n)


rounded = round(in.*10.^n)./10.^n;

end