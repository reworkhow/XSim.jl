# http://weavejl.mpastell.com/stable/
using Weave
cd("/Users/jchen/Dropbox/projects/XSim/tutorials")
filename = "basic"

function render(filename)
    # standalone html
    weave("$filename.jmd")
    # jupyter notebook
    convert_doc("$filename.jmd", "$filename.ipynb")
end


render(filename)