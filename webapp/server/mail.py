"""Outbound email via an authenticated SMTP relay (spec §17: the VM runs no mail
server, so mail leaves through ``QLWEB_SMTP_*``). The transport is injectable so
tests never touch the network -- ``smtp_mailer(cfg)`` returns a plain
``send(to, subject, body)`` callable, and every caller accepts a substitute.

Email addresses are handed to smtplib and NEVER logged here (spec §17: emails
are never logged). This module is import-safe without any SMTP config: the
relay pieces are read from ``cfg`` only when a send actually happens."""
from __future__ import annotations

import smtplib
import ssl
from email.message import EmailMessage

from webapp.server.i18n import t as _t


def smtp_mailer(cfg):
    """Return the default transport: a ``send(to, subject, body)`` callable that
    delivers one plaintext message through the configured STARTTLS relay. The
    connection pieces come from ``cfg`` (host/port/user/pass/from); nothing here
    is evaluated until a message is actually sent, so importing this module
    never requires an SMTP config."""

    def _send(to: str, subject: str, body: str) -> None:
        msg = EmailMessage()
        msg["From"] = cfg.smtp_from
        msg["To"] = to
        msg["Subject"] = subject
        msg.set_content(body)
        with smtplib.SMTP(cfg.smtp_host, cfg.smtp_port, timeout=30) as s:
            s.starttls(context=ssl.create_default_context())
            if cfg.smtp_user:
                s.login(cfg.smtp_user, cfg.smtp_pass)
            s.send_message(msg)

    return _send


def notify_completion(cfg, job, status: str, mailer=None) -> None:
    """Email the big-job requester that their computation finished (or failed),
    linking the permalink page. A no-op when the plaintext address was already
    cleared. The body is rendered in the job's own language; the address is
    passed straight to the transport and never logged."""
    if not job.email:
        return
    send = mailer or smtp_mailer(cfg)
    lang = getattr(job, "lang", "en") or "en"
    url = cfg.public_base_url.rstrip("/") + "/job/" + job.id
    if status == "done":
        subject = _t("mail.done_subject", lang)
        body = _t("mail.done_body", lang).replace("{url}", url)
    else:
        subject = _t("mail.failed_subject", lang)
        body = _t("mail.failed_body", lang).replace("{url}", url)
    send(job.email, subject, body)
