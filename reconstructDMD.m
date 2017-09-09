function [out] = reconstructDMD(modes,freq,b,time)
	out = modes*diag(exp(freq.*time))*b;
	out = real(out);
end
