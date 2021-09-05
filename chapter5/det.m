% função para calcular determinante
classdef calculo_determinante
    methods(Static)
        function [mat_menor]=menor(matriz,i)
            mat_menor=matriz;
            mat_menor(1,:)=[];
            mat_menor(:,i)=[];
        end

        function [det]=calc_det(matriz)
            mat=matriz;
           if size(mat,2)==1
               det=mat(1,1);
               return 
           else
               det=0;
               tam=size(mat,2);
               for x=1:tam
                   det= det + mat(1,x) * ((-1)^(1+x)) * calculo_determinante.calc_det(calculo_determinante.menor(mat,x));
               end
           end
        end
    end
end