function  variance=variance_calcualtion(Array,threshold)
    for i =1:length(threshold)
        if i==1
            subA=Array(Array>threshold(i+1));      
            error = subA - threshold(i);
            mse(i) =  mean(error.^2);
        else if i==length(threshold)
                subA=Array(Array<threshold(i-1));      
                error = subA - threshold(i);
                mse(i) =  mean(error.^2);
            else
                subA=Array((Array<threshold(i-1))&(Array>threshold(i+1)));
                error = subA - threshold(i);
                mse(i) =  mean(error.^2);
            end
        end
    end
    variance=mean(mse);









