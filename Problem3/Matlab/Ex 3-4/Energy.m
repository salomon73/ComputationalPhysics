%% Computing the E-operator %%
function E=E_operator(basis)
    E=cell(1,length(basis));
    for i=1:length(basis)
        E{i}=basis{1,i}*basis{1,i}';
    end
end