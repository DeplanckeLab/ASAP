# ASAP : Automated Single-cell Analysis Pipeline

To access the online tool: <a href="https://asap.epfl.ch">https://asap.epfl.ch</a>.

The web server is composed of 2 main applications:
- [asap_run](https://github.com/DeplanckeLab/asap_run) container, currently in production (v8), which is used by the server to run all the pipelines. It's versionned separately than the asap_web container. It can also be used separately by users to reproduce results outside of the web page.
- [asap_web](https://github.com/fabdavid/asap2_web) container, which is the current web container in production (v8)
- [asap_web](https://github.com/DeplanckeLab/asap_web) container, which is the current web container in development (v9) [Ruby-On-Rails Server code & Docker files](https://github.com/fabdavid/asap2_web)

# Issues
Please post any global issues about the website on this git repo.

Please visit <a href="https://asap.epfl.ch">ASAP</a> for more information.
