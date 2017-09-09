classdef DMD
    properties
        matM1;
        matM2;
        dmdModes;
        dt;
        eigenVal;
        solved;
        b;
        freq;
    end
    methods
        function obj = feedData(obj,fst,scd,dt)
            obj.matM1=fst;
            obj.matM2=scd;
            if(size(obj.matM1)~=size(obj.matM2))
                error('Data matrices must have same size');
            end
            if(~(isequal(obj.matM1(:,2), obj.matM2(:,1)) && ...
                            isequal(obj.matM1(:,end),obj.matM2(:,end-1))))
                error('Data matrices are not consistent');
            end         
            obj.dt=dt;
            obj.solved=false;
        end
        
        function obj = solve(obj)
            [U,Sigma,V] = svd(obj.matM1,'econ');            
            A_bar = (conj(U))' * obj.matM2 * V * pinv(Sigma);
            disp(['Reciprocal Condition No. of A_bar: ', num2str(rcond(A_bar))]);
            [W, eigVal] = eig(A_bar);
            obj.eigenVal = diag(eigVal);
            [obj.eigenVal, si] = sort(obj.eigenVal,'descend');
            W = W(:,si);
%             obj.dmdModes = U*W;
            obj.dmdModes = obj.matM2 * V * pinv(Sigma) * W;
            obj.b = pinv(obj.dmdModes)*obj.matM1(:,1);
            obj.solved = true;            
            obj.freq = log(obj.eigenVal)./obj.dt;
        end
        function [out,error] = reconstruct(obj,smodes,time)
            out = obj.dmdModes(:,smodes)*diag(exp(obj.freq(smodes).*time))*obj.b(smodes);
						out = real(out);
            if(round(time/obj.dt) < size(obj.matM1,2))
                error = (norm(obj.matM2(:,round(time/obj.dt)))-norm(out))/...                    
                norm(obj.matM2(:,round(time/obj.dt)));
            else
                error = 1;
            end
        end
    end
        
end
