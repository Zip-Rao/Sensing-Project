API 参考
========

``sqc`` 的公开 API,按八层 cQED 栈(自底向上)加全局配置模块组织。每个条目展开为一
个页面,列出该子包所导出(其 ``__all__``)的类与函数。

.. note::

   这里只记录 v1 公开面。计划在未来版本提供的能力(以及被排除在某子包 ``__all__``
   之外的名字)是有意省略的;详见路线图。

.. currentmodule:: sqc

.. autosummary::
   :toctree: generated
   :recursive:

   config
   devices
   hardware
   control
   simulation
   experiments
   reconstruction
   calibration
   workflows
