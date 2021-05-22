# http://weavejl.mpastell.com/stable/
using Weave

filename = "simple"

# standalone html
weave("$filename.jmd")

# jupyter notebook
convert_doc("$filename.jmd", "$filename.ipynb")