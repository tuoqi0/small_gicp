# C++ API documentation

See the [README](../README.md) for dependencies, building, native API usage and
tests.

To generate local C++ reference pages, install Doxygen and run:

```sh
cd docs
make cpp
```

The generated [local C++ API](doc_cpp/index.html) is not tracked and is not built
by the normal ROS workspace command. `make mkdocs` builds the documentation site
when MkDocs and its theme are installed.

The [published upstream C++ reference](https://koide3.github.io/small_gicp/doc_cpp/index.html)
and retained [paper](paper.md) describe the upstream publication. Use this fork's
README and installed headers for its native API.
