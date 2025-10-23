function mms = Moments(n,Y_val,mu)

count = 1;
for ell = 0:n
    for k = 1:2*ell+1
        mms(count,:) = mu(ell+1)*Y_val(count,:);
        count = count + 1;
    end
end

end