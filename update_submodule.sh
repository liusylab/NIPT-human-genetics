# add submodule: basevar 
git submodule add https://github.com/ShujiaHuang/BaseVar2.git

# init the submodule and fetch all the files of basevar locally
git submodule update --init --recursive

# update the submodule
git submodule update --remote


# Remove submodule: https://www.geeksforgeeks.org/git/how-to-remove-a-submodule/
git submodule deinit -f BaseVar2  #



rm -rf basevar
