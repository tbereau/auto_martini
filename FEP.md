# Frequently Encountered Problems

## Prediction algorithm

Since ALOGPS, the prediction algorithm for octanol/water partitioning, relies on whole fragments rather than individual atoms, the prediction of certain fragments can pose problem, e.g., small inorganic groups. In this case, `auto_martini` tries to parametrize alternati>
```
; ERROR: no successful mapping found.
; Try running with the '--fpred' and/or '--verbose' options.
```
As mentioned in the error message, an alternative solution consists of relying on an atom-based partitioning coefficient prediction algorithm (Wildman-Crippen), which is less accurate but can predict any fragment.  In case the `--fpred` option is selected, only fragment>

### Boost error

Some versions of Boost will fail to correctly exit at the end of the program, generating such output messages:
```
python: /usr/include/boost/thread/pthread/mutex.hpp:108: boost::mutex::~mutex(): Assertion `!posix::pthread_mutex_destroy(&m)' failed.
[1]    31433 abort (core dumped)  ./auto_martini --smi "N1=C(N)NN=C1N" --mol GUA
```
the results provided by the code are unaffected by this error message. Simply ignore it.

### RDKit outdated

Older RDKit versions will report the following error:
```
[...]
distBdAt = Chem.rdMolTransforms.GetBondLength(conf,i,beadId)
AttributeError: 'module' object has no attribute 'GetBondLength'
```
Simply update your version of RDKit. Most package managers will allow you to do this, unless you've installed RDKit from source.
