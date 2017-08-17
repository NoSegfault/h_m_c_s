

function bal = balance_ncut(deg, f)


    Pf = f - (f'*deg/sum(deg));
    bal = deg'*abs(Pf);
    
end