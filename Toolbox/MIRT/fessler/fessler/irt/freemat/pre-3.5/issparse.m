function y = issparse(x)
y = strcmp(typeof(x), 'sparse');
