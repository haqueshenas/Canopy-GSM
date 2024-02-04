function R = ST3(datavect,m,p)
        vect = datavect(:,1);
        imax = 1;
        for i=1:length(vect)
            if(length(cell2mat(vect(imax)))~=0)
               if(cell2mat(vect(i))>cell2mat(vect(imax)))
                   imax = i;
               end
            else
                imax = i;
            end
        end
      
        
        R = zeros(1,length(vect));
        for i=1:length(vect)
            if(length(cell2mat(vect(i)))~=0)
           if(i<=m)
               Befor = cell2mat(vect(1:i));
               After = cell2mat(vect(i:i+m));
           else
               if(i>length(vect)-m)
                   Befor = cell2mat(vect(i-m:i));
                   After = cell2mat(vect(i:end));
               else
                   Befor = cell2mat(vect(i-m:i));
                   After = cell2mat(vect(i:i+m));
               end
           end
           Befor2 = Befor<1;
           After2 = After<1;
           if((sum(Befor2)==length(Befor2)& (length(Befor2)>m) )...
                   | (sum(After2)==length(After2)& (length(After2)>m) )...
                   ) 
              R(i) = 2; 
           else
               Befor3 = Befor>=1;
               After3 = After>=1;
               if ((length(Befor2)>m  & sum(Befor3)==length(Befor3)) ...
                       | (sum(After3)==length(After3) & ( length(After2)>m))...
                       | (i == imax ))
                   R(i) = 4;
               elseif((i>p & i<imax) ...
                       & sum(Befor3)>=1 ...
                       & sum(After3)>=1 ...
                       & (~(sum(After3)==length(After3)))...
                       & (length(Befor2)>m & length(After2)>m) )
                    R(i) = 3;
               
               elseif(( i > p))
                   R(i) = 5;
               else
                   %% 
                   R(i) = 1;
               end
           end
            else
                R(i)=0;
            end
        end
end